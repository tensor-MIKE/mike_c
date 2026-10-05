#!/usr/bin/env python3
import sys, itertools
from math import ceil, floor, log
import sage.all

class FpEl:
    def __init__(self, n, p, d_word_params, montgomery=True):
        self.n = n
        self.p = p
        self.d_word_params = d_word_params # Dictionnary mapping word sizes to radices
        self.montgomery = montgomery
    def __get_radix(self, word_size, arith=None):
        if arith == "ref" or arith is None:
            if word_size in (32,64):
                return self.d_word_params[word_size]["Radix"]
            else:
                raise ValueError("Word size should be 32 or 64 bits")
        elif arith == "sat64":
            return word_size
        raise ValueError(f'Invalid arithmetic implementation type \"{arith}\"')
    def __get_nlimbs(self, word_size, arith=None):
        if arith == "ref" or arith is None:
            if word_size in (32,64):
                return self.d_word_params[word_size]["Nlimbs"]
            else:
                raise ValueError("Word size should be 32 or 64 bits")
        elif arith == "sat64":
            return ceil(self.d_word_params[word_size]["Nbits"]/word_size)
        raise ValueError(f'Invalid arithmetic implementation type \"{arith}\"')
    def _literal(self, sz, arith=None):
        radix = self.__get_radix(sz, arith=arith)
        l = self.__get_nlimbs(sz, arith=arith)
        # If we're using Montgomery representation, we need to multiply
        # by the Montgomery factor R = 2^lw (l = limb number, w = radix)
        if self.montgomery:
            R = 2**(radix * ceil(log(self.p, 2**radix)))
        else:
            R = 1
        el = (self.n * R) % self.p
        vs = [(int(el) >> radix*i) % 2**radix for i in range(l)]
        return '{' + ', '.join(map(hex, vs)) + '}'
    
class MpEl:
    def __init__(self, x, p):
        self.x = x
        self.p = p
    def _literal(self, sz, arith=None):
        d = sz
        n = (int(self.p).bit_length()//d)+1
        vs =  [(self.x//(2**(d*(i)))%(2**d)) for i in range(n)]
        return '{' + ', '.join(map(hex, vs)) + '}'

class Object:
    def __init__(self, ty, name, obj):
        if(ty=="digit_t[]"):
            self.ty = "digit_t","[]"
        elif '[' in ty:
            idx = ty.index('[')
            depth = ty.count('[]')
            def rec(os, d):
                assert d >= 0
                if not d:
                    return ()
                assert isinstance(os,list) or isinstance(os,tuple)
                r, = {rec(o, d-1) for o in os}
                return (len(os),) + r
            dims = rec(obj, depth)
            self.ty = ty[:idx], ''.join(f'[{d}]' for d in dims)
        else:
            self.ty = ty, ''
        self.name = name
        self.obj = obj

    def _declaration(self):
        return f'extern const {self.ty[0]} {self.name}{self.ty[1]};'

    def _literal(self):
        def rec(obj):
            if isinstance(obj, int):
                if obj < 256: return str(obj)
                else: return hex(obj)
            if isinstance(obj, sage.all.Integer):
                if obj < 256: return str(obj)
                else: return hex(obj)
            if isinstance(obj, FpEl):
                literal = "\n#if 0"
                for sz in (32, 64):
                    literal += f"\n#elif RADIX == {sz}"
                    if sz == 64:
                        literal += "\n#if defined(MIKE_GF_IMPL_SAT64)"
                        literal += f"\n{obj._literal(sz, 'sat64')}"
                        literal += "\n#else"
                        literal += f"\n{obj._literal(sz, 'ref')}"
                        literal += "\n#endif"
                    else:
                        literal += f"\n{obj._literal(sz, 'ref')}"
                return literal + "\n#endif\n"
            if isinstance(obj, MpEl):
                literal = "\n#if 0"
                for sz in (32, 64):
                    literal += f"\n#elif RADIX == {sz}"
                    if sz == 64:
                        literal += "\n#if defined(MIKE_GF_IMPL_SAT64)"
                        literal += f"\n{obj._literal(sz, 'sat64')}"
                        literal += "\n#else"
                        literal += f"\n{obj._literal(sz, 'ref')}"
                        literal += "\n#endif"
                    else:
                        literal += f"\n{obj._literal(sz, 'ref')}"
                return literal + "\n#endif\n"
            if isinstance(obj, list) or isinstance(obj, tuple):
                return '{' + ', '.join(map(rec, obj)) + '}'
            if isinstance(obj, str):
                return obj
            raise NotImplementedError(f'unknown type {type(obj)} in Formatter')
        return rec(self.obj)

    def _definition(self):
        return f'const {self.ty[0]} {self.name}{self.ty[1]} = ' + self._literal() + ';'

class ObjectFormatter:
    def __init__(self, objs):
        self.objs = objs

    def header(self, file=None):
        for obj in self.objs:
            assert isinstance(obj, Object)
            print(obj._declaration(), file=file)

    def implementation(self, file=None):
        for obj in self.objs:
            assert isinstance(obj, Object)
            print(obj._definition(), file=file)
