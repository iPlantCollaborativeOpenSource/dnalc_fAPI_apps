## -> NO USER-SERVICABLE PARTS INSIDE
# Create unique but human comprehensible name for each container

if [ -z "${DOCKER_APP_IMAGE}" ];
then
    echo "Error: DOCKER_APP_IMAGE was not defined"
fi
if [ -z "$HOST_SCRATCH" ];
then
    echo "Warning: HOST_SCRATCH was not defined"
fi

#TTY="-t"
TTY=""
MYUID=$(id -u $USER)
STAMP=$(date +%s)
DOCKER_APP_CONTAINER="app-$STAMP"
HOST_OPTS="-u=$MYUID"
# Launch the app container, persist up to 24h
DOCKER_APP_CREATE="docker run ${HOST_OPTS} -d -v `pwd`:${HOST_SCRATCH}:rw -w ${HOST_SCRATCH} --name ${DOCKER_APP_CONTAINER} ${DOCKER_APP_IMAGE} sleep 86400"

# Export Container ID
export DOCKER_APP_CONTAINER=$( ${DOCKER_APP_CREATE} | cut -c1-12)
# Append this in front of every command
export DOCKER_APP_RUN="docker exec -i ${TTY} $DOCKER_APP_CONTAINER"
## <- NO USER-SERVICABLE PARTS INSIDE
