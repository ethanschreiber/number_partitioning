#!/usr/bin/python

import sys, getopt

def nextChar(c) :
   if c == 'Z' :
      return 'A'
   else :
      return chr(ord(c[-1])+1)
   
def nextString(s) :

   l = list(s)    # Create ist
   idx = len(l)-1

   while idx >= 0 :
      l[idx] = nextChar(l[idx])
      
      if (l[idx] != 'A') :
         break
      
      idx = idx - 1
      
   if (idx == -1) :
      l.insert(0,'A')
   
   return "".join(l)
   
s = 'A'

for i in range(0,54) :
   print s," ",
   s= nextString(s) 
   


