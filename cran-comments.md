---
title: "cran-comments"
output: html_document
---

## R Test environments
* Windows Server 2022 x64, R devel, platform x86_64-w64-mingw32 (64-bit)
* Windows Server 2022 x64, R release, platform x86_64-w64-mingw32
* macOS 12.7.4, R devel
* macOS arm64 14.4.1, R devel
* atlas
* clang-asan
* Linux, R devel

## R CMD check results
- There were no ERRORs, WARNINGs or NOTEs in most of these environments.
- Errors on Linux and Fedora OS when building vignettes (pandoc issue, seems to disappear in CRAN checks according to users)

## Downstream dependencies
There are currently no downstream dependencies for this package.
