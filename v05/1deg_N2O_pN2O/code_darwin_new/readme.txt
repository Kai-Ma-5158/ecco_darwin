
##   添加apN2O输入，准备相关代码  ##  


# =============================在3deg_N2O_20240312/code_darwin_new/基础上进行=============================================

cd /fs2/home/yfeng1982_3/darwin3/1deg_apN2O/code/code_darwin_new/
touch readme.txt

cp /fs2/home/yfeng1982_3/darwin3/3deg_N2O_20240312/code_darwin_new/* .

%%================================ copy the code different from darwin3 ==================================================

cp ../../../pkg/darwin/DARWIN_EXF_PARAMS.h .
cp ../../../pkg/darwin/DARWIN_EXF_FIELDS.h .
cp ../../../pkg/darwin/darwin_exf_load.F .
cp ../../../pkg/darwin/darwin_exf_readparms.F .
cp ../../../pkg/darwin/darwin_exf_init_varia.F .
cp ../../../pkg/darwin/darwin_fields_load.F .
cp ../../../pkg/darwin/darwin_monitor.F .
