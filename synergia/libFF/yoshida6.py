#!/usr/bin/env python

from decimal import Decimal, getcontext

getcontext().prec=30

def e2( n ):
    return Decimal(2)**(Decimal(1)/(Decimal(2)*Decimal(n) + Decimal(1)))

def x0( n ):
    return - e2(n) / (Decimal(2) - e2(n))

def x1( n ):
    return Decimal(1) / (Decimal(2) - e2(n))

c1 = x1(1) * x1(2) / Decimal(2)
c2 = x1(2) * (x0(1) + x1(1)) / Decimal(2)
c3 = x1(1) * (x0(2) + x1(2)) / Decimal(2)
c4 = x0(2) * (x0(1) + x1(1)) / Decimal(2)

d1 = x1(1) * x1(2)
d2 = x0(1) * x1(2)
d3 = x0(2) * x1(1)
d4 = x0(1) * x0(2)

print '    const double c1 = ' + str(c1) + ';'
print '    const double c2 = ' + str(c2) + ';'
print '    const double c3 = ' + str(c3) + ';'
print '    const double c4 = ' + str(c4) + ';'

print '    const double d1 = ' + str(d1) + ';'
print '    const double d2 = ' + str(d2) + ';'
print '    const double d3 = ' + str(d3) + ';'
print '    const double d4 = ' + str(d4) + ';'
