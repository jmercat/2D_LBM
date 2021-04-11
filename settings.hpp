#ifndef SETTINGS_HPP
#define SETTINGS_HPP

#define DYNAMIC_ALLOCATION
#define NAVIERSTOKES
//#define STOKES

#define iterPerCall 50


#define gridSizeX 150
#define gridSizeY 400

#define directionContrast 100

#ifndef EIGEN_STACK_ALLOCATION_LIMIT
// default 131072 == 128 KB, 2097152 = 2MB
#define EIGEN_STACK_ALLOCATION_LIMIT 5097152
#endif

#endif // SETTINGS_HPP
