
chome=$PWD
impl=$PWD/lib/barseqR/barseqRImpl.py
tst=$PWD/test/barseqR_server_test.py
mymod=$PWD/lib/RunDir
util_dir=$PWD/lib/Util
tmp_dir=$PWD/test_local/workdir/tmp/
ui_dir=$PWD/ui/narrative/methods/run_barseqR/
uspec=$PWD/ui/narrative/methods/run_barseqR/spec.json
udisp=$PWD/ui/narrative/methods/run_barseqR/display.yaml
Rdir=$PWD/lib/RunDir/R_dir

#Docker fix
docker run -it -v /var/run/docker.sock:/run/docker.sock alpine chmod g+w /run/docker.sock
