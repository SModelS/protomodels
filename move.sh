#!/bin/sh

TAG="pm200"

git tag -d $TAG
git push origin :$TAG
git tag $TAG
git push origin $TAG
