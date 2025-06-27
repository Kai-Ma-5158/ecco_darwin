

# 在使用mm_100初始文件的基础上，增加apN2O

# ================
# 2. Build executable
# Prerequisite: 1. Get code
cd /fs2/home/yfeng1982_3/darwin3/1deg_apN2O/
mkdir build1 run1
cd build1
rm *


# 最终版(4)
ln -sf ../code/code_darwin_new/* .
rm DIAGNOSTICS_SIZE*
再执行原编译语句
../../tools/genmake2 -mpi -of ../../tools/build_options/linux_amd64_ifort+mpi_ice_eX -mo \ ' ../1deg_code/code_darwin_v4r4_n2o ../1deg_code/code_darwin_v4r4 ../1deg_code/code_v4r4'


make depend
make -j 16

==============
# 3. Instructions for running simulation (1992-2017 period)



（2）-run
cd ../run1
rm -rf *
mkdir -p diags
ln -sf ../build1/mitgcmuv .
cp ../1deg_code/input_v4r4/* .
rm data data.pkg data.diagnostics
cp ../1deg_code/input_darwin_v4r4/* .


INPUTDIR='/fs2/home/yfeng1982_3/Downloads/LLC90/Release4'
ln -s ${INPUTDIR}/input_bin/* .
ln -s ${INPUTDIR}/input_forcing/* .

ln -sf ../apN2O_B2017/* .

ln -sf /fs2/home/yfeng1982_3/Downloads/LLC90/Release5/NOAA_MBL/*1992* .
ln -sf /fs2/home/yfeng1982_3/Downloads/LLC90/Release5/V4r5/* .

mv pickup_ptracers.0000000001.35tra.2.data pickup_ptracers.0000000002.data
mv pickup_ptracers.0000000001.35tra.meta pickup_ptracers.0000000002.meta
rm pickup_ptracers.0000000001.*
cp ../1deg_code/input_darwin_v4r4_new/* .
cp ../code/input_darwin/* .

# qsub job_ECCOV4r4

======提交任务==========
cp ../../1deg/run_llc90*.sh run_llc90_v4r4_n2o.sh
yhbatch run_llc90_v4r4_n2o.sh



=================================== forcing modofy [1992*8 + year+8] ========================






