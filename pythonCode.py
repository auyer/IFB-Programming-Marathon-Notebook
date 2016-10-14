# Fermat Theorem - Teste de Primalidade
def fermat(n):
    if n == 2:
        return True
    if not n & 1:
        return False
    return pow(2, n-1, n) == 1

# Read all input at once
import sys
from itertools import imap # not needed in Python 3, use map instead
data = imap(int, sys.stdin.read().split())
scan = data.next # or data.__next__ if you happen to use Python 3

# print a list of numbers separated by space (or things that can be cast to str)
print " ".join(map(str, L)) # L is the list

# Custom Sort
def qsort(inlist):
    if inlist == []:
        return []
    else:
        pivot = inlist[0]
        lesser = qsort([x for x in inlist[1:] if x < pivot])
        greater = qsort([x for x in inlist[1:] if x >= pivot])
        return lesser + [pivot] + greater

#Save interactiveShell content
import readline
readline.write_history_file('path.py')
