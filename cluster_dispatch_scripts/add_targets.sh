cdp target add --host koustav.pal@vanda.nus.edu.sg \
--remote-root /scratch/projects/CFP03/CFP03-SF-111/koustav.pal/012-Ben-Clarke-ERAP1-KD/ \
--cpus 4 \
--memory 32G \
--time "72:00:00" \
--node 1 \
--scheduler pbs \
--template-file $HOME/hpc_template/pbs_default_template.sh \
vanda_compute

cdp target add --host koustav@patanilab01.lsi.nus.edu.sg \
--remote-root /home/koustav/Projects/012-Ben-Clarke-ERAP1-KD/ \
--cpus 4 \
--memory 8G \
--time "72:00:00" \
--node 1 \
--scheduler sge \
--parallel-environment make \
--template-file $HOME/hpc_template/cluster_dispatch_sge_template.sh \
patani_compute

# cdp target add --host koustav@optimus.bic.nus.edu.sg \
# --remote-root /export2/home/koustav/Projects/012-Ben-Clarke-ERAP1-KD/ \
# --cpus 4 \
# --memory 32G \
# --time "72:00:00" \
# --node 1 \
# --queue all.q \
# --parallel-environment openmpi \
# --scheduler sge \
# --template-file $HOME/hpc_template/cluster_dispatch_sge_template.sh \
# lsi_compute