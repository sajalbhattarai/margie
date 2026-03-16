#!/bin/bash
# Container version of run_with_bvbrc.sh
# The environment is already set up in the Dockerfile (ENV vars)
# So we just exec the command.

# However, we can double check or source if needed. 
# But since we set ENV in Dockerfile, we don't need to source anything.

exec "$@"
