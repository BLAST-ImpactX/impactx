#ifndef FF_MISC_ELEMENTS_H
#define FF_MISC_ELEMENTS_H

#include "ff_element.h"
#include "ff_drift.h"

// =========== monitor =================
class FF_monitor : public FF_drift
{
public:
    FF_monitor() {}
    virtual ~FF_monitor() {}
};

// =========== hmonitor =================
class FF_hmonitor : public FF_drift
{
public:
    FF_hmonitor() { }
    virtual ~FF_hmonitor() {}
};

// =========== vmonitor =================

class FF_vmonitor : public FF_drift
{
public:
    FF_vmonitor() { }
    virtual ~FF_vmonitor() {}
};

// =========== instrument =================

class FF_instrument : public FF_drift
{
public:
    FF_instrument() { }
    virtual ~FF_instrument() {}
};

#endif // FF_MISC_ELEMENTS_H
