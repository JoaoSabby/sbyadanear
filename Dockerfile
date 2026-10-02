# Build the canonical base first:
# docker build -f docker/oraclelinux97-r453/Dockerfile -t sbyadanear-oneapi .
FROM sbyadanear-oneapi
WORKDIR /workspace/sbyadanear
