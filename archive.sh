# git garbage collector reduce .git size
git -c gc.reflogExpire=0 -c gc.reflogExpireUnreachable=0 -c gc.rerereresolved=0 -c gc.rerereunresolved=0 -c gc.pruneExpire=now gc "$@"
tar --exclude='.venv' --exclude='backward_facing_step_compressible/*/*' --exclude='backward_facing_step_compressible/[1-9]*' --exclude='*.tar.gz' --exclude='*.tar' -cvf archive.tar * .git
# extract using tar -xvf archive.tar