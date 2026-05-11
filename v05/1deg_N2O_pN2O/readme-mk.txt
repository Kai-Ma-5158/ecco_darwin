# ========
# 1. Get code
git clone --branch darwin_ckpt68d_at_c66g https://github.com/jahn/darwin3
git clone --depth 1 https://github.com/MITgcm-contrib/ecco_darwin.git
cd darwin3
mkdir build


# ================
# 2. Build executable
cd build
cd /fs2/home/yfeng1982_3/darwin3/1deg_apN2O/
cp ../code_darwin_new/DARWIN_DIAGS.h .
cp ../code_darwin_new/darwin_plankton.F .
cp ../code_darwin_new/darwin_forcing.F .
cp ../code_darwin_new/darwin_diagnos*init* .
cp /fs2/home/yfeng1982_3/darwin3/1deg/1deg_code/code_darwin_v4r4_n2o/darwin_forcing.F .
ln -sf ../code/code_darwin_new/* .
ln -sf ../code/code_darwin_new_add_diags/* .
rm DIAGNOSTICS_SIZE*

../../tools/genmake2 -mpi -of ../../tools/build_options/linux_amd64_ifort+mpi_ice_eX -mo \ ' ../1deg_code/code_darwin_v4r4_n2o ../1deg_code/code_darwin_v4r4 ../1deg_code/code_v4r4'
make depend
make -j 16


# ====================
# 3. Instructions for running simulation (8-year cycled 1992, and 1992-2017 period)
mkdir run
cd ../run
mkdir -p diags diags/monthly diags/daily diags/3hourly diags/budget

#
ln -sf ../build/mitgcmuv .

# input setting
cp ../1deg_code/input_v4r4/* .
rm data data.pkg data.diagnostics
cp ../1deg_code/input_darwin_v4r4/* .
cp ../1deg_code/input_darwin_v4r4_new/* .
cp ../code/input_darwin/* 


INPUTDIR='/fs2/home/yfeng1982_3/Downloads/LLC90/Release4'
# atmospheric forcing
ln -s ${INPUTDIR}/input_forcing/* .
# apN2O forcing
ln -sf ../apN2O_B2017/*.
#apCO2 forcing
ln -sf /fs2/home/yfeng1982_3/Downloads/LLC90/Release5/NOAA_MBL/* .

# init
ln -s ${INPUTDIR}/input_bin/* .
ln -sf /fs2/home/yfeng1982_3/Downloads/LLC90/Release5/V4r5/* .
mv pickup_ptracers.0000000001.35tra.2.data pickup_ptracers.0000000002.data
mv pickup_ptracers.0000000001.35tra.meta pickup_ptracers.0000000002.meta
rm pickup_ptracers.0000000001.*

====== submit ==========
yhbatch run_llc90_v4r4_n2o.sh








