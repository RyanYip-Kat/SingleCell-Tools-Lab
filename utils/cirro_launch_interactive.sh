adata=$1
port=$2
/home/ye/anaconda3/envs/pegasus/bin/cirro launch $adata \
	--host 10.100.44.197 \
	--port $port \
        --no-open	
