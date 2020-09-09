FROM alpine:3.9.2

RUN apk update && apk add --no-cache git bash openssh git-lfs
COPY compose/production/large_data_service/update_data /update_data
RUN chmod +x /update_data