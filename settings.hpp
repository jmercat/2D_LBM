#ifndef SETTINGS_HPP
#define SETTINGS_HPP

#define DYNAMIC_ALLOCATION
#define NAVIERSTOKES
//#define STOKES

#define iterPerCall 20


#define gridSizeX 100
#define gridSizeY 200

#ifndef EIGEN_STACK_ALLOCATION_LIMIT
// default 131072 == 128 KB, 2097152 = 2MB
#define EIGEN_STACK_ALLOCATION_LIMIT 2097152
#endif

#endif // SETTINGS_HPP
