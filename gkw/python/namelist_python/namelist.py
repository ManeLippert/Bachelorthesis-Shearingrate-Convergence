#!/usr/bin/python3

import unittest
try:
    from collections import OrderedDict
except ImportError:
    from .utils import OrderedDict

import re
import sys

def read_namelist_file(filename):
    return Namelist(open(filename, 'r').read())

# trick for py2/3 compatibility
if 'basestring' not in globals():
   basestring = str

class AttributeMapper():
    """
    Simple mapper to access dictionary items as attributes.
    """

    def __init__(self, obj):
        self.__dict__['data'] = obj

    def __getattr__(self, attr):
        if attr in self.data:
            found_attr = self.data[attr]
            if isinstance(found_attr, dict):
                return AttributeMapper(found_attr)
            else:
                return found_attr
        else:
            raise AttributeError

    def __setattr__(self, attr, value):
        if attr in self.data:
            self.data[attr] = value
        else:
            raise NotImplementedError

    def __dir__(self):
        return self.data.keys()

class CaseInsensitiveDict(OrderedDict):
    """
    This is an ordered dictionary which ignores the letter case of the
    keys given (if the respective key is a string).

    Many thanks to user m000 who answered
    https://stackoverflow.com/questions/2082152/case-insensitive-dictionary
    so helpfully.
    """
    
    @classmethod
    def _k(cls, key):
        return key.lower() if isinstance(key, basestring) else key

    def __init__(self, *args, **kwargs):
        super(CaseInsensitiveDict, self).__init__(*args, **kwargs)
        self._convert_keys()
    def __getitem__(self, key):
        return super(CaseInsensitiveDict, self).__getitem__(self.__class__._k(key))
    def __setitem__(self, key, value):
        super(CaseInsensitiveDict, self).__setitem__(self.__class__._k(key), value)
    def __delitem__(self, key):
        return super(CaseInsensitiveDict, self).__delitem__(self.__class__._k(key))
    def __contains__(self, key):
        return super(CaseInsensitiveDict, self).__contains__(self.__class__._k(key))
    def has_key(self, key):
        return super(CaseInsensitiveDict, self).has_key(self.__class__._k(key))
    def pop(self, key, *args, **kwargs):
        return super(CaseInsensitiveDict, self).pop(self.__class__._k(key), *args, **kwargs)
    def get(self, key, *args, **kwargs):
        return super(CaseInsensitiveDict, self).get(self.__class__._k(key), *args, **kwargs)
    def setdefault(self, key, *args, **kwargs):
        return super(CaseInsensitiveDict, self).setdefault(self.__class__._k(key), *args, **kwargs)
    def update(self, E={}, **F):
        super(CaseInsensitiveDict, self).update(self.__class__(E))
        super(CaseInsensitiveDict, self).update(self.__class__(**F))
    def _convert_keys(self):
        for k in list(self.keys()):
            v = super(CaseInsensitiveDict, self).pop(k)
            self.__setitem__(k, v)
    # def to_dict(self):
    #     import copy
    #     ret = {}
    #     for k in self.keys():
    #         if isinstance(self[k], CaseInsensitiveDict):
    #             ret[k] = self[k].to_dict()
    #         else:
    #             ret[k] = copy.deepcopy(self[k])
    #     return ret

class Namelist():
    """
    Parses namelist files in Fortran 90 format.

    Note that while Fortran speaks of "namelists", this module uses
    the term "group" to refer to individual namelists within a file..

    After parsing,

    nlist = namelist_python.read_namelist_file('input.dat')

    recognised groups are accessible through
    the 'groups' attribute (which is a case-insensitive ordered dictionary)

    nlist.groups['mode']['chin']

    or the data attribute

    nlist.data.mode.chin

    """

    def __init__(self, input_str="", parse_strings_unqoted=True):
        self.groups = CaseInsensitiveDict()
        self.parse_strings_unqoted = parse_strings_unqoted

        namelist_start_line_re = re.compile(r'^\s*&(\w+)\s*$')
        # FIXME the end of the namelist does not necessarily have to
        # be in a separate line
        namelist_end_line_re = re.compile(r'^\s*/\s*$')

        # a pattern matching an array of stuff
        a_number = r'[0-9\.\+\-eE]+'

        # a comma-separated list, of elements which either do not
        # contain a comma, or may contain commas inside strings, or
        # may contain commas inside paretheses.
        # At the end of the line, the comma is optional.
        # FIXME strings containing parentheses will cause problems with this expression
        array_re = re.compile(r"(\s*(?:[0-9]+\*)?(?:[^,(\'\"]+|\'[^\']*\'|\"[^\']*\"|[(][^),]+,[^),]+[)])\s*)\s*(?: |,|,?\s*$)")
        self._complex_re = re.compile(r'\s*\([^,]+,[^,]+\)\s*')

        # a pattern to match the non-comment part of a line. This
        # should be able to deal with ! signs inside strings.
        comment_re = re.compile(r"((?:[^\'!]*(?:\'[^\']*\'))*)!.*")

        # match notation for Fortran logicals in namelist files:
        self.logical_true_re = re.compile(r"[^tTfF\']*[tT].*")
        self.logical_false_re = re.compile(r"[^tTfF\']*[fF].*")

        # match abbreviated lists of identical items, like
        # 509*-1.0000000000000000
        # 253*0
        # NEO_EQUIL_PARSE_SP_SEQ=          1,          2, 2*3          , 28*-1
        # 60*"          "
        self.abbrev_list_re = re.compile(r"\s*([0-9]+)\*(.+)\s*")

        # commas at the end of lines seem to be optional
        keyval_line_re = re.compile(r"\s*([\w\(\)]+)\s*=\s*(.*),?")

        # detect index notation for arrays
        array_index_notation_re = re.compile(r"\s*(\w+)\(([0-9])+\)\s*")

        current_group = None
        for line in input_str.split('\n'):
            # remove comments and whitespaces
            line_without_comment = comment_re.sub(r"\1",line).strip()

            if len(line_without_comment) == 0:
                continue

            m = namelist_start_line_re.match(line_without_comment)
            if(m):
                found_group = m.group(1).lower()
                if(current_group is None):
                    if(found_group in self.groups):
                        if(not isinstance(self.groups[found_group],list)):
                            self.groups[found_group] = [self.groups[found_group]]

                        self.groups[found_group].append(CaseInsensitiveDict())
                        current_group = self.groups[found_group][-1]
                        
                    else:
                        self.groups[found_group] = CaseInsensitiveDict()
                        current_group = self.groups[found_group]
                        
                    continue
                else:
                    raise SyntaxError('Namelist %s starts, but namelist %s is not yet complete.' % (found_group,current_group))

            m = namelist_end_line_re.match(line_without_comment)
            if(m):
                if(current_group is not None):
                    current_group = None
                    continue
                else:
                    raise SyntaxError('End of namelist encountered, but there is no corresponding open namelist.')

            # other lines: key = value, or a continuation line
            m = keyval_line_re.match(line_without_comment)
            if(m):
                if(current_group is not None):
                    variable_name = m.group(1)
                    variable_value = m.group(2)

                    # check if this is in array index notation
                    m_ind_notation = array_index_notation_re.match(variable_name)
                    if(m_ind_notation):
                        variable_name = m_ind_notation.group(1)
                        # Fortran indexing is 1-based,
                        # but used Python indexing here
                        index = int(m_ind_notation.group(2))-1
                        #print("index notation: %s %i" % (variable_name, index))
                        if(variable_name not in current_group):
                            current_group[variable_name] = (index+1)*[None]
                        elif(len(current_group[variable_name]) <= index):
                            current_group[variable_name].extend((index+1-len(current_group[variable_name]))*[None])

                    # parse the array with self-crafted regex
                    parsed_list = array_re.findall(variable_value)
                    parsed_list = [self._parse_value(elem) for elem in parsed_list]
                    parsed_list = self._flatten(parsed_list)

                    # if it wasnt for special notations like .false. or 60*'' , one could
                    # simply use a parser from the python standard library for the right
                    # hand side as a whole:
                    #parsed_value = ast.literal_eval(variable_value)
                    try:
                        if(len(parsed_list) == 1):
                            if(m_ind_notation):
                                current_group[variable_name][index] = parsed_list[0]
                            else:
                                current_group[variable_name] = parsed_list[0]
                        else:
                            if(m_ind_notation):
                                current_group[variable_name][index] = parsed_list
                            else:
                                current_group[variable_name] = parsed_list
                    except TypeError:
                        if(m_ind_notation):
                            current_group[variable_name][index] = parsed_list
                        else:
                            current_group[variable_name] = parsed_list
                else:
                    raise SyntaxError('Key %s encountered, but there is no enclosing namelist' % variable_name)
            else:
                warning_text = 'this line could not be parsed: %s' % line_without_comment
                print("WARNING: %s" % warning_text, file=sys.stderr)
                #raise SyntaxError(warning_text)

    def _parse_value(self, variable_value_str):
        """
        Tries to parse a single value, raises a SyntaxError if not successful.
        """
        import ast
        try:
            parsed_value = ast.literal_eval(variable_value_str.strip())

            # use a regex to check if value is a complex number: (1.2 , 3.4)
            # this is needed, because literal_eval parses both "(1.2 , 3.4)"
            # and "1.2 , 3.4" into a tupel with two elements and then one
            # cannot distinguish between a list of two numbers and a single
            # complex number. This makes a difference when it comes to
            # dumping, though.
            if(self._complex_re.match(variable_value_str)):
                parsed_value = complex(parsed_value[0], parsed_value[1])

            try:
                if(isinstance(parsed_value, basestring)):
                    # Fortran strings end with blanks
                    parsed_value = parsed_value.rstrip()
                else:
                    parsed_value = [elem.rstrip() for elem in parsed_value]
            except Exception as err:
                # value is probably just not iterable, or is an iterable of numbers
                pass
        except (ValueError, SyntaxError):

            abbrev_list_match = self.abbrev_list_re.match(variable_value_str)
            if(abbrev_list_match):
                parsed_value = int(abbrev_list_match.group(1)) * [self._parse_value(abbrev_list_match.group(2))]
            elif(self.logical_true_re.match(variable_value_str) and
                 (variable_value_str.lower() in ['true','.true','.true.','t'] or not self.parse_strings_unqoted)):
                parsed_value = True
                if(variable_value_str.lower() not in ['true','.true','.true.','t'] and not self.parse_strings_unqoted):
                    print("WARNING: value %s was parsed to boolean %s" % (variable_value_str, parsed_value), file=sys.stderr)
            elif(self.logical_false_re.match(variable_value_str) and
                 (variable_value_str.lower() in ['false','.false','.false.','f'] or not self.parse_strings_unqoted)):
                parsed_value = False
                if(variable_value_str.lower() not in ['false','.false','.false.','f'] and not self.parse_strings_unqoted):
                    print("WARNING: value %s was parsed to boolean %s" % (variable_value_str, parsed_value), file=sys.stderr)
            else:
                quoted = "'" + variable_value_str.strip()+"'"
                try:
                    parsed_value = ast.literal_eval(quoted)
                    print("WARNING: value %s was treated as %s" % (variable_value_str, quoted), file=sys.stderr)
                except:
                    raise SyntaxError('Right hand side expression could not be parsed. The string is: %s' % (variable_value_str))

                #FIXME distinguish complex scalar and a list of 2 reals
        try:
            if(len(parsed_value) == 1):
                # one gets a list of length 1 if the line ends with a
                # comma, because (4,) for python is a tuple with one
                # element, and (4) is just the scalar 4.
                return parsed_value[0]
            else:
                return parsed_value
        except TypeError:
            return parsed_value
            

    def dump(self, array_inline=True, float_format="%13.5e"):
        lines = []
        for group_name, group_content in self.groups.items():

            group_list = isinstance(group_content,list) and group_content or [group_content]
            for group in group_list:
                lines.append("&%s" % group_name.upper())
                for variable_name, variable_value in group.items():
                    if(isinstance(variable_value, list) or isinstance(variable_value, tuple)):
                        if(array_inline and None not in variable_value):
                            lines.append("%s= %s" % (variable_name, ", ".join([self._format_value(elem, float_format) for elem in variable_value])))
                        else:
                            for n, v in enumerate(variable_value):
                                if(v is not None):
                                    lines.append("%s(%d)= %s" % (variable_name, n+1, self._format_value(v, float_format)))
                    else:
                        lines.append("%s=%s" % (variable_name, self._format_value(variable_value, float_format)))
                lines.append("/")
                lines.append("")

        return "\n".join(lines)

    def _flatten(self, x):
        result = []
        for elem in x:
            if hasattr(elem, "__iter__") and not isinstance(elem, basestring):
                result.extend(self._flatten(elem))
            else:
                result.append(elem)
        return result
    
    def _format_value(self, value, float_format):
        if isinstance(value, bool):
            return value and '.true.' or '.false.'
        elif isinstance(value, int):
            return "%d" % value
        elif isinstance(value, float):
            return float_format % value
        elif isinstance(value, basestring):
            return "'%s'" % value
        elif isinstance(value, complex):
            complex_format = "("+float_format+","+float_format+")"
            return complex_format % (value.real,value.imag)
        else:
            print(value)
            raise Exception("Variable type not understood: type %s" % type(value))

    # create a read-only propery by using property() as a
    # decorator. This function is then the getter function for the
    # .data attribute:
    @property
    def data(self):
        return AttributeMapper(self.groups)

class ParsingTests(unittest.TestCase):

    def __init__(self, methodName='runTest'):
        super().__init__(methodName)
        self.maxDiff = None

    def test_single_value(self):
        input_str = """
        &CCFMSIM_SETUP
        ccfmrad=800.0
        /
        """
        namelist = Namelist(input_str)

        expected_output = {'ccfmsim_setup': { 'ccfmrad': 800. }}

        self.assertEqual(namelist.groups, expected_output)

    def test_multigroup(self):
        input_str = """
        &CCFMSIM_SETUP
        ccfmrad=800.0
        /
        &GROUP2
        R=500.0
        /
        """
        namelist = Namelist(input_str)

        expected_output = {'ccfmsim_setup': { 'ccfmrad': 800. },
                           'group2': { 'r': 500. }}

        self.assertEqual(namelist.groups, expected_output)

    def test_comment(self):
        input_str = """
        ! Interesting comment at the start
        &CCFMSIM_SETUP
        ccfmrad=800.0
        ! And a comment some where in the middle
        /
        &GROUP2
        r=500.0
        /
        """
        namelist = Namelist(input_str)

        expected_output = {'ccfmsim_setup': { 'ccfmrad': 800. },
                           'group2': { 'r': 500. }}

        self.assertEqual(namelist.groups, expected_output)

    def test_array(self):
        input_str = """
        &CCFMSIM_SETUP
        ntrac_picture=4
        var_trac_picture(1)='watcnew'
        des_trac_picture(1)='cloud_water'
        var_trac_picture(2)='watpnew'
        des_trac_picture(2)='rain'
        var_trac_picture(3)='icecnew'
        des_trac_picture(3)='cloud_ice'
        var_trac_picture(4)='granew'
        des_trac_picture(4)='graupel'
        /
        """
        namelist = Namelist(input_str)

        expected_output = {
            'ccfmsim_setup': {
                'ntrac_picture': 4,
                'var_trac_picture': [
                    'watcnew',
                    'watpnew',
                    'icecnew',
                    'granew',
                ],
                'des_trac_picture': [
                    'cloud_water',
                    'rain',
                    'cloud_ice',
                    'graupel',
                ],
            },
        }

        self.assertEqual(dict(namelist.groups), expected_output)

    def test_boolean_sciformat(self):
        input_str = """
        &ATHAM_SETUP

        nz      =300
        zstart  =0.
        ztotal  =15000.
        dzzoom  =50.012345e-15
        kcenter =20
        nztrans =0
        nztrans_boundary =6

        cpumax  =9.e6

        no_uwind=.false.
        no_vwind=.true.
        /
        """
        namelist = Namelist(input_str)

        expected_output = {
            'atham_setup': {
                'nz': 300,
                'zstart': 0.,
                'ztotal': 15000.,
                'dzzoom': 50.012345e-15,
                'kcenter': 20,
                'nztrans': 0,
                'nztrans_boundary': 6,
                'cpumax': 9.e6,
                'no_uwind': False,
                'no_vwind': True,
            }
        }

        self.assertEqual(dict(namelist.groups), expected_output)

    def test_comment_with_forwardslash(self):
        input_str = """
        ! Interesting comment at the start
        &CCFMSIM_SETUP
        ccfmrad=800.0
        ! And a comment some where in the middle/halfway !
        var2=40
        /
        &GROUP2
        r=500.0
        /
        """
        namelist = Namelist(input_str)

        expected_output = {'ccfmsim_setup': { 'ccfmrad': 800., 'var2': 40 },
                           'group2': { 'r': 500. }}

        self.assertEqual(namelist.groups, expected_output)

    def test_inline_array(self):
        input_str = """
        ! can have blank lines and comments in the namelist input file
        ! place these comments between NAMELISTs

        !
        ! not every compiler supports comments within the namelist
        !   in particular vastf90/g77 does not
        !
        ! some will skip NAMELISTs not directly referenced in read
        !&BOGUS rko=1 /
        !
        &TTDATA
        ttreal =  1.,
        ttinteger = 2,
        ttcomplex = (3.,4.),
        ttchar = 'namelist',
        ttbool = T/
        &AADATA
        aareal =  1.  1.  2.  3.,
        aainteger = 2 2 3 4,
        aacomplex = (3.,4.) (3.,4.) (5.,6.) (7.,7.),
        aachar = 'namelist' 'namelist' 'array' ' the lot',
        aabool = T T F F/
        &XXDATA
        xxreal =  1.,
        xxinteger = 2,
        xxcomplex = (3.,4.)/! can have blank lines and comments in the namelist input file
        """

        expected_output = {
            'ttdata': {
                'ttreal': 1.,
                'ttinteger': 2,
                'ttcomplex': 3. + 4.j,
                'ttchar': 'namelist',
                'ttbool': true,
            },
            'aadata': {
                'aareal': [1., 1., 2., 3.,],
                'aainteger': [2, 2, 3, 4],
                'aacomplex': [3.+4.j, 3.+4.j, 5.+6.j, 7.+7.j],
                'aachar': ['namelist', 'namelist', 'array', ' the lot'],
                'aabool': [True, True, False, False],
            },
            'xxdata': {
                'xxreal': 1.,
                'xxinteger': 2.,
                'xxcomplex': 3.+4.j,
            },
        }

        namelist = Namelist(input_str)

        self.assertEqual(dict(namelist.groups), expected_output)

    def test_single_value(self):
        input_str = """&CCFMSIM_SETUP
ccfmrad=  8.00000e+02
/
"""
        namelist = Namelist(input_str)

        self.assertEqual(namelist.dump(), input_str)

    def test_multigroup(self):
        input_str = """&CCFMSIM_SETUP
ccfmrad=  8.00000e+02
/

&GROUP2
r=  5.00000e+02
/
"""
        namelist = Namelist(input_str)

        self.assertEqual(namelist.dump(), input_str)


    def test_array(self):
        input_str = """&CCFMSIM_SETUP
var_trac_picture(1)= 'watcnew'
var_trac_picture(2)= 'watpnew'
var_trac_picture(3)= 'icecnew'
var_trac_picture(4)= 'granew'
des_trac_picture(1)= 'cloud_water'
des_trac_picture(2)= 'rain'
des_trac_picture(3)= 'cloud_ice'
des_trac_picture(4)= 'graupel'
/
"""
        namelist = Namelist(input_str)

        self.assertEqual(namelist.dump(array_inline=False), input_str)

    def test_inline_array(self):
        input_str = """&AADATA
aacomplex= (3.000000,4.000000) (3.000000,4.000000) (5.000000,6.000000) (7.000000,7.000000)
bbcomplex= (  3.00000e+00,  4.00000e+00), (  3.00000e+00,  4.00000e+00), (  5.00000e+00,  6.00000e+00), (  7.00000e+00,  7.00000e+00)
/
"""

        namelist = Namelist(input_str)
        expected_output = {
            'aadata': {
                'aacomplex': [complex(3.,4.),complex(3.,4.),complex(5.,6.),complex(7.,7.),],
                'bbcomplex': [complex(3.,4.),complex(3.,4.),complex(5.,6.),complex(7.,7.),]
            }}
        self.assertEqual(dict(namelist.groups), expected_output)

    def test_inline_array_dump(self):
         input_str = """&DATA
temp_coef=    1.00000000,   5.00000000,   0.18077700,   0.28930000,   0.01450000,
dens_coef=2.00000000,3.00000000,4.18077700,5.28930000,6.01450000
/
"""
         namelist = Namelist(input_str)
         expected_output = {
            'data': {
                'temp_coef': [1.00000000,   5.00000000,   0.18077700,   0.28930000,   0.01450000,],
                'dens_coef': [2.00000000,3.00000000,4.18077700,5.28930000,6.01450000]
            }}
         self.assertEqual(dict(namelist.groups), expected_output)

         expected_output = """&DATA
temp_coef=   1.00000e+00,   5.00000e+00,   1.80777e-01,   2.89300e-01,   1.45000e-02
dens_coef=   2.00000e+00,   3.00000e+00,   4.18078e+00,   5.28930e+00,   6.01450e+00
/
"""

         self.assertEqual(namelist.dump(), expected_output)


if __name__=='__main__':
    unittest.main()
