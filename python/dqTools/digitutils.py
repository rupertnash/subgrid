"""Some digit fiddling tools used to make Sun Grid Engine array
jobs more workable.

In retrospect I would have produced a look up table:
   int ID -> parameters

n is an integer eg. n = 123 = abc (a=1, b=2, c=3)

split(abc) -> [a, b, c]

combine([a, b, c]) -> abc

"""
def split(n, base=10):
    quot, rem = divmod(n, base)
    if quot == 0:
        return [rem]
    else:
        return split(quot,base) + [rem]

def combine(list, base=10):
    if len(list):
        return list[-1] + base*combine(list[:-1], base)
    else:
        return 0

