# for dry-run add -n as an option
rsync -ahmr \
    --verbose \
    --info=progress2 \
    --exclude=archive \
    --exclude=archive2 \
    --exclude='*/restart/' \
    --include='*/' \
    --include-from=include.rsync  \
    --exclude='*' \
    --exclude=jobs \
    -e 'ssh -p 61022' \
    jnattila@gateway.flatironinstitute.org:ceph/runko-gpu/projects/qed-turb/ .
    #jnattila@rusty:ceph/runko/projects/shocks3d/ .


    #jnattila@beskow.pdc.kth.se:/cfs/klemming/nobackup/j/jnattila/runko/projects/turbulence/ .


