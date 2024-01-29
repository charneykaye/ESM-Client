# Use an existing lightweight Linux image as a base
FROM alpine:latest

# Install curl
RUN apk --no-cache add curl

# Set up a default command to run when the container starts
# This can be overridden with a different command when starting the container
CMD ["curl", "--help"]