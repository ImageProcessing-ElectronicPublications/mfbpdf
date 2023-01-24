# MFBpdf

**MFBpdf** is a simple project for easy converting PNM to (MASK+FG+BG)-pdf.
It uses [libtiff](https://github.com/vadz/libtiff) and [libjpeg](https://github.com/LuaDist/libjpeg) for all technichal work and compression.
The breakdown of the image into components is done using [DjVuL](https://github.com/plzombie/depress/issues/2) and [DjVuL wiki](https://sourceforge.net/p/imthreshold/wiki/DjVuL/?version=3).

See [MFBpdf demo](https://github.com/ImageProcessing-ElectronicPublications/mfbpdf-demo).

## Install

### install library

* [libtiff](https://github.com/vadz/libtiff)
* [libjpeg](https://github.com/LuaDist/libjpeg)

### build

Type:

```shell
$ make
```

## Usage

```shell
mfbpdf [options] input.ppm output.tif output.pdf
```

where options are:

```
 -a #      anisotropic regulator DjVuL {0.0}
 -b #      BG and FG downsample {3}
 -c #      contrast regulator DjVuL {0.0}
 -d #      DPI pdf and tiff {300}
 -f #      FG more downsample {2}
 -l #      level DjVuL (0 - auto) {0}
 -o #      overlay blocks DjVuL {0.5}
 -q #      jpeg quality {75}
 -r        rewrite tiff
 -x #      linear regulator DjVuL {*1.0}
 -y #      linear regulator DjVuL {+0.0}
 -z #      black mode
```
You can use [Netpbm](https://sourceforge.net/projects/netpbm/) or any other similar tool to obtain `PNM` from other format. Tiff image for recognition. For example, [tesseract](https://github.com/tesseract-ocr/tesseract)

## DjVuL description.

The [base of the algorithm](https://sourceforge.net/p/imthreshold/wiki/DjVuL/?version=3) was obtained in 2016 by studying the works [monday2000](http://djvu-soft.narod.ru/) and adapting them to Linux.
The prerequisite was the [BookScanLib](http://djvu-soft.narod.ru/bookscanlib/) project  and the algorithm [DjVu Thresholding Binarization](http://djvu-soft.narod.ru/bookscanlib/034.htm).
This algorithm embodied good ideas, but had a recursive structure, was a "function with discontinuities" and had a hard color limit.
The result of this algorithm, due to the indicated shortcomings and the absence of regulators, was doubtful.
After careful study, all the foundations of the specified algorithm were rejected.
The new algorithm is based on levels instead of recursion, a smooth weight function is used instead of a "discontinuous" one, no color restriction, BG/FG selection controls are enabled.
The new algorithm allowed not only to obtain a much more adequate result, but also gave derivative functions: image division into BG/FG according to the existing mask.

## Links

* [pdfbeads](https://github.com/ifad/pdfbeads)
* [qpdf](https://github.com/qpdf/qpdf)
* [pdfcook](https://github.com/ksharindam/pdfcook)
* [tesseract](https://github.com/tesseract-ocr/tesseract)
* [hocr-tools](https://github.com/ocropus/hocr-tools)
* [imthreshold](https://github.com/ImageProcessing-ElectronicPublications/imthreshold)
* [aithreshold](https://github.com/ImageProcessing-ElectronicPublications/aithreshold)
