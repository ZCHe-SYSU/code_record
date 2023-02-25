#download from TCGA database
read -p "Please enter working directory		" x
echo "get into working directory"
cd /Users/hexintong/Documents/WD
mkdir -p run
cd run 
if [ ! -e $x ];then  
	mkdir $x
	echo "creating directory……"
fi

cd $x
if [ ! -e GDCdata ];then
	echo "downloading data......"
	Rscript /Users/hexintong/Documents/WD/script/rscript/download.R $x
fi

if [ ! -e Clinical/vital_status ];then
	echo "obtaining OS&PFS......"
	mkdir -p Clinical
	cd /Users/hexintong/Documents/WD/run/LUAD/GDCdata/TCGA-LUAD/harmonized/Clinical/Clinical_Supplement
	ls */*.xml|grep clinical|awk -F "/" '{print $2}' > name
	rm -rf data
	mkdir data
	for i in `cat name`;do cat */$i |awk -F ">" '{print $2}'|awk -F "</" '{print $2"\t"$1}'|sort|uniq |awk -F ":" '{print $2}'|grep -v '^$' > data/$i ;done
	cd data
	grep days_to_new_tumor_event_after_initial_treatment *|awk -F "." '{print $3"\t"$4}'|cut -f 1,3|sort -nrk2 > ../../../../../../Clinical/new
	grep days_to_death *|awk -F "." '{print $3"\t"$4}'|cut -f 1,3|sort -nrk2 > ../../../../../../Clinical/days_to_death
	grep days_to_last_followup *|awk -F "." '{print $3"\t"$4}'|cut -f 1,3|sort -nrk2 > ../../../../../../Clinical/days_to_last_followup
	grep pathologic_stage  *|awk -F "." '{print $3"\t"$4}'|cut -f 1,3|sort -nrk2 > ../../../../../../Clinical/pathologic_stage
	grep gender *|awk -F "." '{print $3"\t"$4}'|cut -f 1,3|sort -nrk2 > ../../../../../../Clinical/gender
	grep vital_status *|awk -F "." '{print $3"\t"$4}'|cut -f 1,3|sort -nrk2 > ../../../../../../Clinical/vital_status
	cd /Users/hexintong/Documents/WD/run/$x
fi

if [ ! -e dds.Rdata ];then
	echo "data preprocessing......"
	Rscript /Users/hexintong/Documents/WD/script/rscript/preprocessing.R 
fi

if [ ! -e Homo_sapiens.GRCh38.101.abinitio.gtf.gz ];then
	echo "downloading gene notations......"
	wget ftp://ftp.ensembl.org/pub/release-101/gtf/homo_sapiens/Homo_sapiens.GRCh38.101.abinitio.gtf.gz
fi

if [ ! -e dds_cox.csv ];then
	echo "COX univariate analysis......"
	Rscript /Users/hexintong/Documents/WD/script/rscript/cox_survival_analyse.R
fi

if [ ! -e dds_cox.csv-1 ];then
	echo "randomforest……"
	Rscript /Users/hexintong/Documents/WD/script/rscript/randomforest.R $x
fi

if [ ! -e dds_cox.csv-1 ];then
	echo "GO analysising......"
	Rscript /Users/hexintong/Documents/WD/script/rscript/restomerge.R
	python  /Users/hexintong/Documents/WD/script/pscript/merge.py
	Rscript /Users/hexintong/Documents/WD/script/rscript/goanalyse.R
fi
