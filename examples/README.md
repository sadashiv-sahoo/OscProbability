
## Standard Oscillation Examples

## How to Use :
1. cd v1.0.0/
2. cd PREM/
3. chmod +x PREM.sh
4. ./PREM
   PREM484.dat (user's profile)
   [after successful completion]

(optional steps, one may skip to step 10)

5. cd ../
6. pwd $(copy the dir path)
7. edit ~/.bashrc
8. export Oscv1 =$(copied the dir path) (save it !!)
9. source ~/.bashrc
10. cd example/Test$
11. cmake -DSET_TARGET=Test$ $Oscv1        (if step 5 to 9 followed)
    cmake -DSET_TARGET=Test$ -DROOTSYS=ON $Oscv1 (if one uses ROOT6)
    cmake -DSET_TARGET=Test$ ../../v1.0.0   (if step 5 to 9 skipped)

(Assuming all 9 steps are followed)
## Test1
  Vacuum Oscillation !!!
  cmake -DSET_TARGET=Test1 $Oscv1
  make -jN
  ./Test1

## Test2
  Constant density Matter Oscillation !!!
  cmake -DSET_TARGET=Test2 $Oscv1
  make -jN
  ./Test2

## Test3
  Profiled density(DUNE) Matter Oscillation !!!
  cmake -DSET_TARGET=Test3 $Oscv1
  make -jN
  ./Test3

## Test4
  Earth Profiled density (484 Layers) Matter Oscillation !!!
  cmake -DSET_TARGET=Test4 $Oscv1
  make -jN
  ./Test4




