#ifndef DEBUG_HPP
#define DEBUG_HPP

#ifdef DEBUG
    #define ONLYDEBUG(p) p
    #define ONLYRELEASE(p)
#else
    #define ONLYDEBUG(p)
    #define ONLYRELEASE(p) p
#endif // DEBUG_HPP

#endif
