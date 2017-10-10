"""
memory_manager.py - generators for defining / allocating / transfering memory
for kernel creation
"""

from __future__ import division

from .. import utils

from string import Template
import numpy as np
import loopy as lp
import re
import yaml
from enum import Enum
import six
from ..core import array_creator as arc
import re


class memory_type(Enum):
    m_constant = 0,
    m_local = 1,
    m_global = 1

    def __int__(self):
        return self.value


class memory_limits(object):
    """
    Helps determine whether a kernel is using too much constant / shared memory,
    etc.

    Properties
    ----------
    arrays: dict
            A mapping of :class:`loopy.TemporaryVariable` or :class:`loopy.GlobalArg`
            to :class:`memory_type`, representing the types of all arrays to be
            included
    limits: dict
        A dictionary with keys 'shared', 'constant' and 'local' indicated the
        total amount of each memory type available on the device
    """

    def __init__(self, lang, arrays, limits):
        """
        Initializes a :class:`memory_limits`
        """
        self.lang = lang
        self.arrays = arrays
        self.limits = limits

    def can_fit(self, type=memory_type.m_constant, with_type_changes={}):
        """
        Determines whether the supplied :param:`arrays` of type :param:`type`
        can fit on the device

        Parameters
        ----------
        type: :class:`memory_type`
            The type of memory to use
        with_type_changes: dict
            Updates to apply to :prop:`arrays`

        Returns
        -------
        can_fit: int
            The maximum number of times these arrays can be fit in memory.
            - For global memory, this determines the number of initial conditions
            that can be evaluated.
            - For shared and constant memory, this determines whether this data
            can be fit
        """

        if self.lang == 'c':
            return True

        # filter arrays by type
        arrays = self.arrays[type]
        arrays = [a for a in arrays if not any(
            a in v for k, v in six.iteritems(with_type_changes) if k != type)]

        per_ic = 0
        static = 0
        for array in arrays:
            size = 1
            is_ic_dep = False
            for s in array.shape:
                if arc.problem_size.name in str(s):
                    # mark as dependent on # of initial conditions
                    is_ic_dep = True
                    # get the floor div (if any)
                    floor_div = re.search('// (\d+)', str(s))
                    if floor_div:
                        floor_div = int(floor_div.group(1))
                        size /= floor_div
                    continue
                # update size
                size *= s
            # and convert size
            size *= array.dtype.itemsize

            # update counter
            if is_ic_dep:
                per_ic += size
            else:
                static += size

        # if no per_ic, just divide by 1
        if per_ic == 0:
            per_ic = 1

        # return the number of times we can fit these array
        return np.maximum(np.floor((self.limits[type] - static) / per_ic), 0)

    @staticmethod
    def get_limits(loopy_opts, arrays, input_file=""):
        """
        Utility method to load shared / constant memory limits from a file or
        :mod:`pyopencl` as needed

        Parameters
        ----------
        loopy_opts: :class:`loopy_options`
            If loopy_opts.lang == 'opencl', pyopencl will be used to fill in any
            missing limits
        arrays: dict
            A mapping of :class:`memory_type` to :class:`loopy.TemporaryVariable` or
            :class:`loopy.GlobalArg`, representing the types of all arrays to be
            included
        input_file: str
            The path to a yaml file with specified limits (keys should include
            'local' and 'constant', and 'global')

        Returns
        -------
        limits: :class:`memory_limits`
            An initialized :class:`memory_limits` that can determine the total
            'global', 'constant' and 'local' memory available on the device
        """
        limits = {}
        if loopy_opts.lang == 'opencl':
            try:
                limits.update({
                    memory_type.m_global: loopy_opts.device.global_mem_size,
                    memory_type.m_constant:
                        loopy_opts.device.max_constant_buffer_size,
                    memory_type.m_local: loopy_opts.device.local_mem_size})
            except AttributeError:
                pass
            if input_file:
                with open(input_file, 'r') as file:
                    # load from file
                    lims = yaml.load(file.read())
                    mtype = utils.EnumType(memory_type)
                    limits = {}
                    for key, value in six.iteritems(lims):
                        # check in memory type
                        if not key.lower() in [mt.name.lower()[2:]
                                               for mt in memory_type]:
                            msg = ', '.join([t.name.lower() for t in memory_type])
                            msg = '{0}: use one of {1}'.format(memory_type.name, msg)
                            raise Exception(msg)
                        key += 'm_'
                        # update with enum
                        limits[mtype(key)] = value

        return memory_limits(loopy_opts.lang, arrays,
                             {k: v for k, v in six.iteritems(limits)})


class memory_manager(object):

    """
    Aids in defining & allocating arrays for various languages
    """

    def __init__(self, lang):
        """
        Parameters
        ----------
        lang : str
            The language used in this memory initializer
        """
        self.arrays = []
        self.in_arrays = []
        self.out_arrays = []
        self.host_constants = []
        self.lang = lang
        self.memory_types = {'c': 'double*',
                             'opencl': 'cl_mem'}
        self.type_map = {np.dtype('int32'): 'int',
                         np.dtype('float64'): 'double'}
        self.host_langs = {'opencl': 'c',
                           'c': 'c'}
        self.alloc_templates = {'opencl': Template(
            '${name} = clCreateBuffer(context, ${memflag}, '
            '${buff_size}, NULL, &return_code)'),
                                'c': Template(
            '${name} = (double*)malloc(${buff_size})')}
        self.copy_in_templates = {'opencl': Template(
            'clEnqueueWriteBuffer(queue, ${name}, CL_TRUE, 0, '
            '${buff_size}, ${host_buff}, 0, NULL, NULL)'),
                                  'c': Template(
            'memcpy(${name}, ${host_buff}, ${buff_size})')}
        self.copy_out_templates = {'opencl': Template(
            'clEnqueueReadBuffer(queue, ${name}, CL_TRUE, 0, '
            '${buff_size}, ${host_buff}, 0, NULL, NULL)'),
                                   'c': Template(
            'memcpy(${host_buff}, ${name}, ${buff_size})')}
        self.host_constant_template = Template(
            'const ${type} h_${name} [${size}] = {${init}}'
            )
        self.memset_templates = {'opencl': Template(
            """
            #if CL_LEVEL >= 120
                clEnqueueFillBuffer(queue, ${name}, ${fill_value}, ${fill_size}, 0,
                    ${buff_size}, 0, NULL, NULL)
            #else
                clEnqueueWriteBuffer(queue, ${name}, CL_TRUE, 0, ${buff_size},
                    zero, 0, NULL, NULL)
            #endif
            """
            ),
            'c': Template('memset(${name}, 0, ${buff_size})')
        }
        self.free_template = {'opencl': Template('clReleaseMemObject(${name})'),
                              'c': Template('free(${name})')}

    def get_check_err_call(self, call, lang=None):
        if lang is None:
            lang = self.lang
        if lang == 'opencl':
            return Template('check_err(${call})').safe_substitute(call=call)
        else:
            return call

    def add_arrays(self, arrays=[], in_arrays=[], out_arrays=[]):
        """
        Adds arrays to the manager

        Parameters
        ----------
        arrays : list of :class:`lp.GlobalArg`
            The arrays to declare
        in_arrays : list of str
            The array names that form the input to this kernel
        out_arrays : list of str
            The array names that form the output of this kernel

        Returns
        -------
        None
        """
        self.arrays.extend(
            [x for x in arrays if not isinstance(x, lp.ValueArg)])
        self.in_arrays.extend(in_arrays)
        self.out_arrays.extend(out_arrays)

    @property
    def host_arrays(self):
        def _set_sort(arr):
            return sorted(set(arr), key=lambda x: arr.index(x))
        return _set_sort(self.in_arrays + self.out_arrays)

    @property
    def host_lang(self):
        return self.host_langs[self.lang]

    def get_defns(self):
        """
        Returns the definition strings for this memory manager's arrays

        Parameters
        ----------
        None

        Returns
        -------
        defn_str : str
            A string of global memory definitions
        """

        def __add(arraylist, lang, prefix, defn_list):
            for arr in arraylist:
                defn_list.append(self.memory_types[lang] + ' ' + prefix +
                                 arr + utils.line_end[lang])

        defns = []
        # get all 'device' defns
        __add([x.name for x in self.arrays], self.lang, 'd_', defns)

        # return defn string
        return '\n'.join(sorted(set(defns)))

    def get_host_constants(self):
        """
        Returns allocations of initialized constant variables on the host.
        These result when we run out of __constant memory on the device, and must
        migrate to passing constant __global args.

        Parameters
        ----------
        None

        Returns
        -------
        alloc_str : str
            The string of memory allocations
        """

        def _handle_type(arr):
            return lp.types.to_loopy_type(arr.dtype).numpy_dtype

        def _stringify(arr):
            return ', '.join(['{}'.format(x) for x in arr.initializer])

        return '\n'.join([self.host_constant_template.safe_substitute(
            name=x.name,
            type=self.type_map[_handle_type(x)],
            size=str(np.prod(x.shape)),
            init=_stringify(x))
            + utils.line_end[self.lang] for x in self.host_constants])

    def get_mem_allocs(self, alloc_locals=False):
        """
        Returns the allocation strings for this memory manager's arrays

        Parameters
        ----------
        alloc_locals : bool
            If true, only define host arrays

        Returns
        -------
        alloc_str : str
            The string of memory allocations
        """

        def __get_alloc_and_memset(dev_arr, lang, prefix):
            # if we're allocating inputs, call them something different
            # than the input args from python
            post_fix = '_local' if alloc_locals else ''

            in_host_const = any(dev_arr.name == y.name for y in self.host_constants)

            if lang == self.host_lang and in_host_const:
                return ''

            # get name
            name = prefix + dev_arr.name + post_fix
            # if it's opencl, we need to declare the buffer type
            memflag = None
            if lang == 'opencl':
                memflag = 'CL_MEM_READ_WRITE' if not any(
                    x.name == dev_arr.name for x in self.host_constants
                    ) else 'CL_MEM_READ_ONLY'

            # generate allocs
            alloc = self.alloc_templates[lang].safe_substitute(
                name=name,
                memflag=memflag,
                buff_size=self._get_size(dev_arr))

            if alloc_locals:
                # add a type
                alloc = self.memory_types[self.host_lang] + ' ' + alloc

            # generate allocs
            return_list = [alloc]

            # add error checking
            if lang == 'opencl':
                return_list.append(self.get_check_err_call('return_code'))

            if not in_host_const:
                # add the memset
                return_list.append(
                    self.get_check_err_call(
                        self.memset_templates[lang].safe_substitute(
                            name=name,
                            buff_size=self._get_size(dev_arr),
                            fill_value='&zero',  # fill for OpenCL kernels
                            fill_size='sizeof(double)',  # fill type
                            ), lang=lang))

            # return
            return '\n'.join([r + utils.line_end[lang] for r in return_list] +
                             ['\n'])

        to_alloc = [
            next(x for x in self.arrays if x.name == y) for y in self.host_arrays
            ] if alloc_locals else self.arrays
        prefix = 'h_' if alloc_locals else 'd_'
        lang = self.lang if not alloc_locals else self.host_lang
        alloc_list = [__get_alloc_and_memset(
            arr, lang, prefix) for arr in to_alloc]

        return '\n'.join(alloc_list)

    def _get_size(self, arr, subs_n='problem_size'):
        size = arr.shape
        nsize = []
        str_size = []
        skip = []
        # remove 'problem_size' from shape if present, as it's baked into the
        # various defns
        for i, x in enumerate(size):
            s = str(x)
            # check for non-integer sizes
            if 'problem_size' in s:
                str_size.append(subs_n)
                if s != 'problem_size':
                    # it's a floor division thing, need to do some cleanup here
                    vsize = re.search(
                        r'\(-1\)\*\(\(\(-1\)\*problem_size\) // (\d+)\)', s)
                    # make sure we found the vector width
                    assert vsize is not None and len(vsize.groups()) > 0, (
                        'Unknown size for array {}, :{}'.format(arr.name, s))
                    vsize = int(vsize.group(1))
                    assert vsize in size, (
                        'Unknown vector-width: {}, found for array {}'.format(
                            vsize, arr.name))
                    # and eliminate it from the list
                    skip.append(size.index(vsize))
                # skip this entry
                skip.append(i)

        # find all numerical sizes
        nsize = [size[i] for i in range(len(size)) if i not in skip]
        if nsize:
            nsize = str(np.cumprod(nsize, dtype=np.int32)[-1])
            str_size.append(nsize)
        # multiply by the item size
        str_size.append(str(arr.dtype.itemsize))
        return ' * '.join(str_size)

    def _mem_transfers(self, to_device=True, host_constants=False):
        if not host_constants:
            # get arrays
            arr_list = self.in_arrays if to_device else self.out_arrays
            # filter out host constants
            arr_list = [x for x in arr_list if not any(
                y.name == x for y in self.host_constants)]
            arr_maps = {x: next(y for y in self.arrays if x == y.name)
                        for x in arr_list}
        else:
            arr_list = [x.name for x in self.host_constants]
            arr_maps = {x.name: x for x in self.host_constants}

        templates = self.copy_in_templates if to_device else self.copy_out_templates

        return '\n'.join([
            self.get_check_err_call(templates[self.lang].safe_substitute(
                name='d_' + arr, host_buff='h_' + arr,
                buff_size=self._get_size(arr_maps[arr]))) +
            utils.line_end[self.lang] for arr in arr_list])

    def get_mem_transfers_in(self):
        """
        Generates the memory transfers into the device before kernel execution

        Parameters
        ----------
        None

        Returns
        -------
        mem_transfer_in : str
            The string to perform the memory transfers before execution
        """

        return self._mem_transfers(to_device=True)

    def get_mem_transfers_out(self):
        """
        Generates the memory transfers into the back to the host after kernel
        execution

        Parameters
        ----------
        None

        Returns
        -------
        mem_transfer_out : str
            The string to perform the memory transfers back to the host after
            execution
        """

        return self._mem_transfers(to_device=False)

    def get_host_constants_in(self):
        """
        Generates the memory transfers of the host constants
        into the device before kernel execution

        Parameters
        ----------
        None

        Returns
        -------
        mem_transfer_in : str
            The string to perform the memory transfers before execution
        """

        return self._mem_transfers(to_device=True, host_constants=True)

    def get_mem_frees(self, free_locals=False):
        """
        Returns code to free the allocated buffers

        Parameters
        ----------
        free_locals : bool
            If true, we're freeing the local versions of the host arrays

        Returns
        -------
        mem_free_str : str
            The generated code
        """

        if not free_locals:
            # device memory
            frees = [self.get_check_err_call(
                self.free_template[self.lang].safe_substitute(name='d_' + arr.name))
                     for arr in self.arrays]
        else:
            frees = [self.free_template[self.host_lang].safe_substitute(
                name='h_' + arr + '_local') for arr in self.host_arrays if not any(
                        x.name == arr for x in self.host_constants)]

        return '\n'.join([x + utils.line_end[self.lang] for x in sorted(set(frees))])
