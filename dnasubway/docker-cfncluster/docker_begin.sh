## -> NO USER-SERVICABLE PARTS INSIDE
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
UUID=$(uuidgen -r)
HOSTNAME=$(uname -n)
HOST_OPTS="-u=$MYUID"
# Create unique but human comprehensible name for each container
DOCKER_APP_CONTAINER="appc-${UUID}"
PERSIST=86400

# Logging
echo "hostname: $HOSTNAME" >> docker.log
echo "username: $USER" >> docker.log
echo "userid: $MYUID" >> docker.log
echo "container_name: $DOCKER_APP_CONTAINER" >> docker.log

# Refresh the image
docker pull ${DOCKER_APP_IMAGE} | grep "Status" >> docker.log

# Launch the app container, persist up to 24h
DOCKER_APP_CREATE="docker run ${HOST_OPTS} -d -v `pwd`:${HOST_SCRATCH}:rw -w ${HOST_SCRATCH} --name ${DOCKER_APP_CONTAINER} ${DOCKER_APP_IMAGE} sleep ${PERSIST}"

# Export Container ID
export DOCKER_APP_CONTAINER=$( ${DOCKER_APP_CREATE} | cut -c1-12)
echo "container_id: $DOCKER_APP_CONTAINER" >> docker.log
# Append this in front of every command
export DOCKER_APP_RUN="docker exec -i ${TTY} $DOCKER_APP_CONTAINER"
## <- NO USER-SERVICABLE PARTS INSIDE
