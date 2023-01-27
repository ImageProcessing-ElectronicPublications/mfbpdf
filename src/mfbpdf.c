/* $Id$ */

/*
 * Copyright (c) 1991-1997 Sam Leffler
 * Copyright (c) 1991-1997 Silicon Graphics, Inc.
 *
 * Permission to use, copy, modify, distribute, and sell this software and
 * its documentation for any purpose is hereby granted without fee, provided
 * that (i) the above copyright notices and this permission notice appear in
 * all copies of the software and related documentation, and (ii) the names of
 * Sam Leffler and Silicon Graphics may not be used in any advertising or
 * publicity relating to the software without the specific, prior written
 * permission of Sam Leffler and Silicon Graphics.
 *
 * THE SOFTWARE IS PROVIDED "AS-IS" AND WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY
 * WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.
 *
 * IN NO EVENT SHALL SAM LEFFLER OR SILICON GRAPHICS BE LIABLE FOR
 * ANY SPECIAL, INCIDENTAL, INDIRECT OR CONSEQUENTIAL DAMAGES OF ANY KIND,
 * OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
 * WHETHER OR NOT ADVISED OF THE POSSIBILITY OF DAMAGE, AND ON ANY THEORY OF
 * LIABILITY, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE
 * OF THIS SOFTWARE.
 */

#include <time.h>
#include <tiffio.h>
#include <jpeglib.h>

#include "mfbpdf.h"
#define DJVUL_IMPLEMENTATION
#define THRESHOLD_IMPLEMENTATION
#include "djvul.h"
#include "threshold.h"

static void usage(void)
{
    char* stuff[] =
    {
        "usage: mfbpdf [options] input.ppm input|output.tif output.pdf",
        "where options are:",
        " -a #      anisotropic regulator DjVuL {0.0}",
        " -b #      BG and FG downsample {3}",
        " -c #      contrast regulator DjVuL {0.0}",
        " -d #      DPI pdf and tiff {300}",
        " -f #      FG more downsample {2}",
        " -l #      level DjVuL (0 - auto) {0}",
        " -o #      overlay blocks DjVuL {0.5}",
        " -q #      jpeg quality {75}",
        " -r        rewrite tiff",
        " -s #      sensitivity for sauvola and blur {0.2}",
        " -t #      threshold: djvul, bimod, sauvola, blur {djvul}",
        " -x #      linear regulator DjVuL {*1.0}",
        " -y #      linear regulator DjVuL {+0.0}",
        " -z        black mode",
        "",
        NULL
    };

    char buf[BUFSIZ];
    int i;

    setbuf(stderr, buf);
    fprintf(stderr, "%s\n\n", TIFFGetVersion());
    fprintf(stderr, "%s version %s\n\n", "MFBpdf", MFBPDF_VERSION);
    for (i = 0; stuff[i] != NULL; i++)
        fprintf(stderr, "%s\n", stuff[i]);
    exit(-1);
}

static int thresholdoptions(char* opt)
{
	int threshold = TDJVUL;
	
    if (streq(opt, "bimod"))
        threshold = TBIMOD;
    else if (streq(opt, "sauvola"))
        threshold = TSAUVOLA;
    else if (streq(opt, "blur"))
        threshold = TBLUR;

    return (threshold);
}

static void pnmbad(char* file)
{
    fprintf(stderr, "%s: Not a PPM file.\n", file);
    exit(-2);
}

static tmsize_t multiply_ms(tmsize_t m1, tmsize_t m2)
{
    tmsize_t bytes = m1 * m2;

    if (m1 && bytes / m1 != m2)
        bytes = 0;

    return bytes;
}

static unsigned long int jencode(unsigned char **jpeg, unsigned char *buf, int width, int height, int pixelFormat, int quality, int jpegcs, int progressive, int optimize, int subsample)
{
    unsigned long int jpegSize = 0;
    struct jpeg_compress_struct cinfo;
    struct jpeg_error_mgr jerr;
    JSAMPROW row_pointer[1];
    int row_stride = width * ((pixelFormat == JCS_RGB) ? 3 : 1);

    cinfo.err = jpeg_std_error(&jerr);

    jpeg_create_compress(&cinfo);

    // Set destination
    jpeg_mem_dest(&cinfo, jpeg, &jpegSize);

    // Set options
    cinfo.image_width = width;
    cinfo.image_height = height;
    cinfo.input_components = (pixelFormat == JCS_RGB) ? 3 : 1;
    cinfo.in_color_space = pixelFormat;

    jpeg_set_defaults(&cinfo);

    if (optimize)
        cinfo.optimize_coding = TRUE;

    if (optimize && !progressive)
    {
        cinfo.scan_info = NULL;
        cinfo.num_scans = 0;
    }

    if (!optimize && progressive)
    {
        jpeg_simple_progression(&cinfo);
    }

    jpeg_set_quality(&cinfo, quality, TRUE);
    jpeg_set_colorspace (&cinfo, jpegcs);

    if (subsample == 1)
    {
        cinfo.comp_info[0].h_samp_factor = 1;
        cinfo.comp_info[0].v_samp_factor = 1;
        cinfo.comp_info[1].h_samp_factor = 1;
        cinfo.comp_info[1].v_samp_factor = 1;
        cinfo.comp_info[2].h_samp_factor = 1;
        cinfo.comp_info[2].v_samp_factor = 1;
    }

    // Start the compression
    jpeg_start_compress(&cinfo, TRUE);

    // Process scanlines one by one
    while (cinfo.next_scanline < cinfo.image_height)
    {
        row_pointer[0] = &buf[cinfo.next_scanline * row_stride];
        (void) jpeg_write_scanlines(&cinfo, row_pointer, 1);
    }

    jpeg_finish_compress(&cinfo);
    jpeg_destroy_compress(&cinfo);

    return jpegSize;
}

int main(int argc, char* argv[])
{
    unsigned char *buf = NULL;
    unsigned char *bufpbm = NULL;
    unsigned char *bufpnm = NULL;
    bool *bufmask = NULL;
    unsigned char *bufbg = NULL;
    unsigned char *buffg = NULL;
    unsigned char *jcompressed = NULL;
    unsigned long int jsize = 0, tifflen;
    tsize_t kd, ki;
    tsize_t linebytes = 0;
    tsize_t pbmlinebytes = 0;
    tsize_t datasize = 0;
    tsize_t rowsize = 0;
    tsize_t stripsize = 0;
    tsize_t stripcount = 0;
    tsize_t bufferoffset = 0;
    tsize_t read = 0;
    uint16 compression = COMPRESSION_CCITTFAX4;
    uint16 photometric = PHOTOMETRIC_MINISWHITE;
    uint16 spp = 1;
    uint16 bpp = 8;
    uint16 sppm = 1;
    uint16 bppm = 1;
    uint16 bwmw = 0;
    uint16 bwmb = 1;
    TIFF *outtiff;
    FILE *fin, *fout, *ftmp;
    int maskexist = 0;
    unsigned int bgs = 3;
    unsigned int fgs = 2;
    unsigned int level = 0;
    int wbmode = 1;
    int tmode = TDJVUL;
    float doverlay = 0.5f;
    float anisotropic = 0.0f;
    float contrast = 0.0f;
    float fbscale = 1.0f;
    float delta = 0.0f;
    float sensitivity = 0.2f;
    int retiff = 0;
    unsigned int h, w, hm, wm, wd, hbg, wbg,  hfg, wfg;
    unsigned int sc, prec, y, x, d, k;
    int jpegcs = JCS_YCbCr;
    int jpegpf = JCS_RGB;
    int jquality = 75;
    int jsubsample = 0;
    int jprog = 0;
    int jopt = 1;
    int dpi = 300;
    int hpi, wpi, lhi, lwi, phi, pwi, xnum;
    float hp, wp;
    unsigned long int xref[13];
    struct tm* currenttime;
    time_t timenow;
    char *infile;
    int c;
#if !HAVE_DECL_OPTARG
    extern int optind;
    extern char* optarg;
#endif

    if (argc < 4)
    {
        fprintf(stderr, "%s: Too few arguments\n", argv[0]);
        usage();
    }
    while ((c = getopt(argc, argv, "a:b:c:d:f:l:o:q:rs:t:x:y:z")) != -1)
        switch (c)
        {
        case 'a':
            anisotropic = atof(optarg);
            break;
        case 'b':
            bgs = atoi(optarg);
            break;
        case 'c':
            contrast = atof(optarg);
            break;
        case 'd':
            dpi = atoi(optarg);
            break;
        case 'f':
            fgs = atoi(optarg);
            break;
        case 'l':
            level = atoi(optarg);
            break;
        case 'o':
            doverlay = atof(optarg);
            break;
        case 'q':
            jquality = atoi(optarg);
            break;
        case 'r':
            retiff = !retiff;
            break;
        case 's':
            sensitivity = atof(optarg);
            break;
        case 't':
			tmode = thresholdoptions(optarg);
            break;
        case 'x':
            fbscale = atof(optarg);
            break;
        case 'y':
            delta = atof(optarg);
            break;
        case 'z':
            wbmode = !wbmode;
            break;
        case 'h':
        case '?':
            usage();
        }

    if (optind + 3 < argc)
    {
        fprintf(stderr, "%s: Too many arguments\n", argv[0]);
        usage();
    }

    /*
     * If only one file is specified, read input from
     * stdin; otherwise usage is: ppm2tiff input output.
     */
    if (argc - optind > 1)
    {
        infile = argv[optind++];
        fin = fopen(infile, "rb");
        if (fin == NULL)
        {
            fprintf(stderr, "%s: Can not open.\n", infile);
            return (-1);
        }
    }
    else
    {
        infile = "<stdin>";
        fin = stdin;
#if defined(HAVE_SETMODE) && defined(O_BINARY)
        setmode(fileno(stdin), O_BINARY);
#endif
    }

    if (fgetc(fin) != 'P')
        pnmbad(infile);
    switch (fgetc(fin))
    {
    case '4':           /* it's a PBM file */
        bpp = 1;
        spp = 1;
        jpegpf = JCS_GRAYSCALE;
        jpegcs = JCS_GRAYSCALE;
        break;
    case '5':           /* it's a PGM file */
        bpp = 8;
        spp = 1;
        jpegpf = JCS_GRAYSCALE;
        jpegcs = JCS_GRAYSCALE;
        break;
    case '6':           /* it's a PPM file */
        bpp = 8;
        spp = 3;
        jpegpf = JCS_RGB;
        break;
    default:
        pnmbad(infile);
    }

    /* Parse header */
    while(1)
    {
        if (feof(fin))
            pnmbad(infile);
        c = fgetc(fin);
        /* Skip whitespaces (blanks, TABs, CRs, LFs) */
        if (strchr(" \t\r\n", c))
            continue;

        /* Check for comment line */
        if (c == '#')
        {
            do
            {
                c = fgetc(fin);
            }
            while(!(strchr("\r\n", c) || feof(fin)));
            continue;
        }

        ungetc(c, fin);
        break;
    }

    switch (bpp)
    {
    case 1:
        if (fscanf(fin, " %u %u", &w, &h) != 2)
            pnmbad(infile);
        if (fgetc(fin) != '\n')
            pnmbad(infile);
        break;
    case 8:
        if (fscanf(fin, " %u %u %u", &w, &h, &prec) != 3)
            pnmbad(infile);
        if (fgetc(fin) != '\n' || prec != 255)
            pnmbad(infile);
        break;
    }
    switch (bpp)
    {
    case 1:
        /* if round-up overflows, result will be zero, OK */
        linebytes = (multiply_ms(spp, w) + (8 - 1)) / 8;
        break;
    case 8:
        linebytes = multiply_ms(spp, w);
        break;
    }
    pbmlinebytes = (multiply_ms(1, w) + (8 - 1)) / 8;
    bufpbm = (unsigned char *)_TIFFmalloc(pbmlinebytes);
    if (bufpbm == NULL)
    {
        fprintf(stderr, "%s: Not enough memory\n", infile);
        exit(-2);
    }
    bufpnm = (unsigned char *)_TIFFmalloc(w * h * spp);
    if (bufpnm == NULL)
    {
        fprintf(stderr, "%s: Not enough memory\n", infile);
        exit(-2);
    }
    bufmask = (bool *)_TIFFmalloc(w * h);
    if (bufmask == NULL)
    {
        fprintf(stderr, "%s: Not enough memory\n", infile);
        exit(-2);
    }
    hbg = (h + bgs - 1) / bgs;
    wbg = (w + bgs - 1) / bgs;
    hfg = (hbg + fgs - 1) / fgs;
    wfg = (wbg + fgs - 1) / fgs;
    bufbg = (unsigned char *)_TIFFmalloc(wbg * hbg * spp);
    if (bufbg == NULL)
    {
        fprintf(stderr, "%s: Not enough memory\n", infile);
        exit(-2);
    }
    buffg = (unsigned char *)_TIFFmalloc(wbg * hbg * spp);
    if (buffg == NULL)
    {
        fprintf(stderr, "%s: Not enough memory\n", infile);
        exit(-2);
    }
    buf = (unsigned char *)_TIFFmalloc(linebytes);
    if (buf == NULL)
    {
        fprintf(stderr, "%s: Not enough memory\n", infile);
        exit(-2);
    }
    kd = 0;
    for (y = 0; y < h; y++)
    {
        if (fread(buf, linebytes, 1, fin) != 1)
        {
            fprintf(stderr, "%s: scanline %lu: Read error.\n", infile, (unsigned long) y);
            break;
        }
        if (bpp == 1)
        {
            for (x = 0; x < w; x++)
            {
                k = x >> 3;
                sc = ((buf[k] & (0x80 >> (x & 0x07))) == 0) ? 0 : 1;
                bufpnm[kd] = (sc == 1) ? 0 : 255;
                bufmask[kd] = sc;
                kd++;
            }
        }
        else
        {
            ki = 0;
            for (x = 0; x < w; x++)
            {
                for (d = 0; d < spp; d++)
                {
                    bufpnm[kd] = buf[ki];
                    kd++;
                    ki++;
                }
            }
        }
    }
    if (fin != NULL)
        fclose(fin);

    ftmp = fopen(argv[optind], "r");
    if (ftmp != NULL)
    {
        maskexist = 1;
        fclose(ftmp);
    }
    if (maskexist > 0)
    {
        outtiff = TIFFOpen(argv[optind], "r");
        if (outtiff == NULL)
        {
            fprintf(stderr, "Can't open input file %s for reading\n", argv[optind]);
            exit(-2);
        }
        TIFFGetField(outtiff, TIFFTAG_IMAGEWIDTH, &wm);
        TIFFGetField(outtiff, TIFFTAG_IMAGELENGTH, &hm);
        if ((hm != h) || (wm != w))
        {
            fprintf(stderr, "Bad size mask %s: %ux%u\n", argv[optind], wm, hm);
            (void) TIFFClose(outtiff);
            exit(-2);
        }

        if(TIFFGetField(outtiff, TIFFTAG_COMPRESSION, &compression) == 0)
        {
            fprintf(stderr, "No support for %s with compression type %u: not configured\n", argv[optind], compression);
            (void) TIFFClose(outtiff);
            exit(-2);
        }
        if( TIFFIsCODECConfigured(compression) == 0)
        {
            fprintf(stderr, "No support for %s with compression type %u: not configured\n", argv[optind], compression);
            (void) TIFFClose(outtiff);
            exit(-2);
        }
        TIFFGetFieldDefaulted(outtiff, TIFFTAG_BITSPERSAMPLE, &bppm);
        switch(bppm)
        {
        case 0:
            fprintf(stderr, "Image %s has 0 bits per sample, assuming 1\n", argv[optind]);
            bppm=1;
            break;
        case 1:
            break;
        default:
            fprintf(stderr, "Image %s not mask\n", argv[optind]);
            (void) TIFFClose(outtiff);
            exit(-2);
        }
        TIFFGetFieldDefaulted(outtiff, TIFFTAG_SAMPLESPERPIXEL, &sppm);
        switch(sppm)
        {
        case 1:
            break;
        default:
            fprintf(stderr, "Image %s not mask\n", argv[optind]);
            (void) TIFFClose(outtiff);
            exit(-2);
        }
        if(TIFFGetField(outtiff, TIFFTAG_PHOTOMETRIC, &photometric) == 0)
        {
            fprintf(stderr, "No support for %s with no photometric interpretation tag\n", argv[optind]);
            (void) TIFFClose(outtiff);
            exit(-2);
        }
        switch(photometric)
        {
        case PHOTOMETRIC_MINISWHITE:
            break;
        case PHOTOMETRIC_MINISBLACK:
            bwmw = 1;
            bwmb = 0;
            break;
        default:
            fprintf(stderr, "Image %s not mask\n", argv[optind]);
            (void) TIFFClose(outtiff);
            exit(-2);
            break;
        }
        datasize = TIFFTileSize(outtiff);
        unsigned char *buffer = NULL;
        buffer = (unsigned char*)_TIFFmalloc(datasize);
        if (buffer == NULL)
        {
            fprintf(stderr, "Can't allocate %lu bytes of memory for image, %s\n", datasize, argv[optind]);
            (void) TIFFClose(outtiff);
            exit(-2);
        }
        memset(buffer, 0, datasize);
        rowsize = TIFFTileRowSize(outtiff);
        stripsize = TIFFStripSize(outtiff);
        stripcount = TIFFNumberOfStrips(outtiff);
        for(y = 0; y < stripcount; y++)
        {
            read = TIFFReadEncodedStrip(outtiff, y, (tdata_t) &buffer[bufferoffset], stripsize);
            if(read == -1)
            {
                fprintf(stderr, "Error on decoding strip %u of %s\n", y, argv[optind]);
                (void) TIFFClose(outtiff);
                exit(-2);
            }
            bufferoffset += read;
        }
        kd = 0;
        ki = 0;
        for (y = 0; y < h; y++)
        {
            for (x = 0; x < w; x++)
            {
                k = x >> 3;
                sc = ((buffer[ki + k] & (0x80 >> (x & 0x07))) == 0) ? bwmw : bwmb;
                bufmask[kd] = sc;
                kd++;
            }
            ki += rowsize;
        }
        maskexist = 1;
        if (buffer)
            _TIFFfree(buffer);
        if(outtiff != NULL)
            TIFFClose(outtiff);
    }

    if (bpp > 1)
    {
        if (maskexist == 0)
        {
            if(!(level = ImageDjvulSelect(bufpnm, bufmask, bufbg, buffg, w, h, spp, bgs, level, wbmode, doverlay, anisotropic, contrast, fbscale, delta, (float)(bgs * fgs), sensitivity, tmode)))
            {
                fprintf(stderr, "ERROR: not complite DjVuL\n");
                exit(-3);
            }
        }
        else
        {
            if(!(level = ImageDjvulGround(bufpnm, bufmask, bufbg, buffg, w, h, spp, bgs, level, doverlay)))
            {
                fprintf(stderr, "ERROR: not complite DjVuL ground\n");
                exit(-3);
            }
        }
        if (fgs > 1)
        {
            fgs = ImageFGdownsample(buffg, wbg, hbg, spp, fgs);
        }
    }

    if ((maskexist == 0) || (retiff > 0))
    {
        outtiff = TIFFOpen(argv[optind], "w");
        if (outtiff == NULL)
        {
            fprintf(stderr, "Not read %s\n", argv[optind]);
            exit(-2);
        }
        TIFFSetField(outtiff, TIFFTAG_IMAGEWIDTH, (uint32) w);
        TIFFSetField(outtiff, TIFFTAG_IMAGELENGTH, (uint32) h);
        TIFFSetField(outtiff, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
        TIFFSetField(outtiff, TIFFTAG_SAMPLESPERPIXEL, 1);
        TIFFSetField(outtiff, TIFFTAG_BITSPERSAMPLE, 1);
        TIFFSetField(outtiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
        TIFFSetField(outtiff, TIFFTAG_PHOTOMETRIC, photometric);
        TIFFSetField(outtiff, TIFFTAG_COMPRESSION, compression);
        TIFFSetField(outtiff, TIFFTAG_ROWSPERSTRIP, h);
        if (linebytes == 0)
        {
            fprintf(stderr, "%s: scanline size overflow\n", infile);
            (void) TIFFClose(outtiff);
            exit(-2);
        }
        if (dpi > 0)
        {
            TIFFSetField(outtiff, TIFFTAG_XRESOLUTION, dpi);
            TIFFSetField(outtiff, TIFFTAG_YRESOLUTION, dpi);
            TIFFSetField(outtiff, TIFFTAG_RESOLUTIONUNIT, RESUNIT_INCH);
        }

        for (y = 0; y < h; y++)
        {
            for (x = 0; x < w; x++)
            {
                k = x >> 3;
                bufpbm[k] = 0;
            }
        }
        kd = 0;
        for (y = 0; y < h; y++)
        {
            wd = (8 - w % 8) % 8;
            for (x = 0; x < w; x++)
            {
                k = x >> 3;
                kd++;
                bufpbm[k] <<= 1;
                bufpbm[k] |= ((bufmask[kd]) ? 1 : 0);
            }
            k = (w - 1) >> 3;
            bufpbm[k] <<= wd;
            if (TIFFWriteScanline(outtiff, bufpbm, y, 0) < 0)
                break;
        }
        (void) TIFFClose(outtiff);
        if (buf)
            _TIFFfree(buf);
    }
    if (bufmask)
        _TIFFfree(bufmask);
    if (bufpnm)
        _TIFFfree(bufpnm);

    fin = fopen(argv[optind],"rb");
    if (fin == NULL)
    {
        fprintf(stderr, "%s: Can not open.\n", argv[optind]);
        return (-1);
    }
    fseek(fin, 0l, SEEK_END);
    tifflen = ftell(fin) - 8;
    rewind(fin);
    fseek(fin, 8l, SEEK_SET);
    buf = (unsigned char *)_TIFFmalloc(tifflen + 1);
    if (buf == NULL)
    {
        fprintf(stderr, "%s: Not enough memory\n", argv[optind]);
        (void) TIFFClose(outtiff);
        exit(-2);
    }
    if (fread(buf, tifflen, 1, fin) != 1)
    {
        fprintf(stderr, "Read error tiff.\n");
        exit(-3);
    }
    if (fin != NULL)
        fclose(fin);

    optind++;
    fout = fopen(argv[optind], "wb");
    if (fout == NULL)
    {
        fprintf(stderr, "%s: Can not open.\n", argv[optind]);
        return (-1);
    }
    hp = (float)h * 72.0f / (float)dpi;
    wp = (float)w * 72.0f / (float)dpi;
    hpi = (int)(hp * 100.0f + 0.5f);
    wpi = (int)(wp * 100.0f + 0.5f);
    lhi = phi = 1;
    while (hpi > phi)
    {
        phi *= 10;
        lhi++;
    }
    lhi--;
    lwi = pwi = 1;
    while (wpi > pwi)
    {
        pwi *= 10;
        lwi++;
    }
    lwi--;
    xnum = 0;
    xref[xnum] = ftell(fout);
    xnum++;
    fwrite("%PDF-1.4\n%\342\343\317\323\n", 15, 1, fout);
    xref[xnum] = ftell(fout);
    xnum++;
    fwrite("1 0 obj\n", 8, 1, fout);
    fwrite("\n<<\n/Kids [2 0 R]\n/Type /Pages\n/Count 1\n>>\nendobj\n", 50, 1, fout);
    fwrite("\n", 1, 1, fout);
    xref[xnum] = ftell(fout);
    xnum++;
    fwrite("2 0 obj\n", 8, 1, fout);
    fwrite("\n<<\n/pdftk_PageNum 1\n/Contents 3 0 R\n/Type /Page\n/Resources \n<<\n/ProcSet [/PDF /ImageC]\n/ExtGState 4 0 R\n/XObject 5 0 R\n>>\n/Parent 1 0 R\n/MediaBox [0 0 ", 152, 1, fout);
    fprintf(fout, "%.2f %.2f", wp, hp);
    fwrite("]\n>>\nendobj\n", 12, 1, fout);
    fwrite("\n", 1, 1, fout);
    xref[xnum] = ftell(fout);
    xnum++;
    fwrite("3 0 obj\n", 8, 1, fout);
    fwrite("\n<<\n/Length ", 12, 1, fout);
    if (bpp > 1)
    {
        fprintf(fout, "%d", (33 + lhi + 3 + (lwi + 1 + lhi) + 3 + lwi + 13 + lwi + 5 + lhi + 29 + lhi + 3 + (lwi + 1 + lhi) + 3 + lwi + 13 + lwi + 5 + lhi + 30));
    }
    else
    {
        fprintf(fout, "%d", (33 + lhi + 3 + (lwi + 1 + lhi) + 3 + lwi + 13 + lwi + 5 + lhi + 30));
    }
    fwrite("\n>>\nstream\n", 11, 1, fout);
    fwrite("q 0.01 0 0 0.01 0 0 cm\nq\n1 0 m\n1 ", 33, 1, fout);
    fprintf(fout, "%d", hpi);
    fwrite(" l\n", 3, 1, fout);
    fprintf(fout, "%d %d", wpi, hpi);
    fwrite(" l\n", 3, 1, fout);
    fprintf(fout, "%d", wpi);
    fwrite(" 0 l\nh\nW n\nq ", 13, 1, fout);
    fprintf(fout, "%d", wpi);
    fwrite(" 0 0 ", 5, 1, fout);
    fprintf(fout, "%d", hpi);
    if (bpp > 1)
    {
        fwrite(" 0 0 cm\n/R7 Do\nQ\nQ\nq\n1 0 m\n1 ", 29, 1, fout);
        fprintf(fout, "%d", hpi);
        fwrite(" l\n", 3, 1, fout);
        fprintf(fout, "%d %d", wpi, hpi);
        fwrite(" l\n", 3, 1, fout);
        fprintf(fout, "%d", wpi);
        fwrite(" 0 l\nh\nW n\nq ", 13, 1, fout);
        fprintf(fout, "%d", wpi);
        fwrite(" 0 0 ", 5, 1, fout);
        fprintf(fout, "%d", hpi);
        fwrite(" 0 0 cm\n/R9 Do\nQ\nQ\n/R10 gs\nQ\n\n", 30, 1, fout);
    }
    else
    {
        fwrite(" 0 0 cm\n/R8 Do\nQ\nQ\n/R10 gs\nQ\n\n", 30, 1, fout);
    }
    fwrite("endstream\nendobj\n", 17, 1, fout);
    fwrite("\n", 1, 1, fout);
    xref[xnum] = ftell(fout);
    xnum++;
    fwrite("4 0 obj\n", 8, 1, fout);
    fwrite("\n<<\n/R10 6 0 R\n>>\nendobj\n", 25, 1, fout);
    fwrite("\n", 1, 1, fout);
    xref[xnum] = ftell(fout);
    xnum++;
    fwrite("5 0 obj\n", 8, 1, fout);
    if (bpp > 1)
    {
        fwrite("\n<<\n/R7 7 0 R\n/R8 8 0 R\n/R9 9 0 R\n>>\nendobj\n", 44, 1, fout);
    }
    else
    {
        fwrite("\n<<\n/R8 8 0 R\n>>\nendobj\n", 24, 1, fout);
    }
    fwrite("\n", 1, 1, fout);
    xref[xnum] = ftell(fout);
    xnum++;
    fwrite("6 0 obj\n", 8, 1, fout);
    fwrite("\n<<\n/Type /ExtGState\n/TR /Identity\n>>\nendobj\n", 45, 1, fout);
    fwrite("\n", 1, 1, fout);
    if (bpp > 1)
    {
        xref[xnum] = ftell(fout);
        xnum++;
        jsize = jencode(&jcompressed, buffg, wfg, hfg, jpegpf, jquality, jpegcs, jprog, jopt, jsubsample);
        fwrite("7 0 obj\n", 8, 1, fout);
        fwrite("\n<<\n/ColorSpace ", 16, 1, fout);
        if (spp == 1)
        {
            fwrite("/DeviceGray", 11, 1, fout);
        }
        else
        {
            fwrite("/DeviceRGB", 10, 1, fout);
        }
        fwrite("\n/Subtype /Image\n/Height ", 25, 1, fout);
        fprintf(fout, "%d", hfg);
        fwrite("\n/Filter /DCTDecode\n/Width ", 27, 1, fout);
        fprintf(fout, "%d", wfg);
        fwrite("\n/BitsPerComponent 8\n/Length ", 29, 1, fout);
        fprintf(fout, "%ld", jsize);
        fwrite("\n>>\nstream\n", 11, 1, fout);
        fwrite(jcompressed, jsize, 1, fout);
        fwrite("\nendstream\nendobj\n", 18, 1, fout);
        fwrite("\n", 1, 1, fout);
    }
    xref[xnum] = ftell(fout);
    xnum++;
    fwrite("8 0 obj\n", 8, 1, fout);
    fwrite("\n<<\n/ColorSpace /DeviceGray\n/Subtype /Image\n/Height ", 52, 1, fout);
    fprintf(fout, "%d", h);
    fwrite("\n/Filter /CCITTFaxDecode\n/DecodeParms\n<<\n/Columns ", 50, 1, fout);
    fprintf(fout, "%d", w);
    fwrite("\n/K -1\n>>\n/Width ", 17, 1, fout);
    fprintf(fout, "%d", w);
    fwrite("\n/BitsPerComponent 1\n/Length ", 29, 1, fout);
    fprintf(fout, "%ld", tifflen);
    fwrite("\n>>\nstream\n", 11, 1, fout);
    fwrite(buf, tifflen, 1, fout);
    fwrite("\nendstream\nendobj\n", 18, 1, fout);
    fwrite("\n", 1, 1, fout);
    if (bpp > 1)
    {
        xref[xnum] = ftell(fout);
        xnum++;
        jsize = jencode(&jcompressed, bufbg, wbg, hbg, jpegpf, jquality, jpegcs, jprog, jopt, jsubsample);
        fwrite("9 0 obj\n", 8, 1, fout);
        fwrite("\n<<\n/ColorSpace ", 16, 1, fout);
        if (spp == 1)
        {
            fwrite("/DeviceGray", 11, 1, fout);
        }
        else
        {
            fwrite("/DeviceRGB", 10, 1, fout);
        }
        fwrite("\n/Subtype /Image\n/Height ", 25, 1, fout);
        fprintf(fout, "%d", hbg);
        fwrite("\n/Filter /DCTDecode\n/Width ", 27, 1, fout);
        fprintf(fout, "%d", wbg);
        fwrite("\n/SMask 8 0 R\n/BitsPerComponent 8\n/Length ", 42, 1, fout);
        fprintf(fout, "%ld", jsize);
        fwrite("\n>>\nstream\n", 11, 1, fout);
        fwrite(jcompressed, jsize, 1, fout);
        fwrite("\nendstream\nendobj\n", 18, 1, fout);
        fwrite("\n", 1, 1, fout);
    }
    xref[xnum] = ftell(fout);
    xnum++;
    fwrite("10 0 obj\n", 9, 1, fout);
    fwrite("\n<<\n/Type /Catalog\n/Pages 1 0 R\n>>\nendobj\n", 42, 1, fout);
    fwrite("\n", 1, 1, fout);
    xref[xnum] = ftell(fout);
    xnum++;
    fwrite("11 0 obj\n", 9, 1, fout);
    fwrite("\n<<\n/Creator (mfbpdf 1.0)\n", 26, 1, fout);
    fwrite("/Producer (MFB template)\n", 25, 1, fout);
    fwrite("/CreationDate ", 14, 1, fout);
    if (time(&timenow) == (time_t) -1)
    {
        timenow = (time_t) 0;
    }
    currenttime = gmtime(&timenow);
    fprintf(fout, "(D:%.4d%.2d%.2d%.2d%.2d%.2d)\n",
            (currenttime->tm_year + 1900) % 65536,
            (currenttime->tm_mon + 1) % 256,
            (currenttime->tm_mday) % 256,
            (currenttime->tm_hour) % 256,
            (currenttime->tm_min) % 256,
            (currenttime->tm_sec) % 256);
    fwrite("/ModDate ", 9, 1, fout);
    fprintf(fout, "(D:%.4d%.2d%.2d%.2d%.2d%.2d)\n",
            (currenttime->tm_year + 1900) % 65536,
            (currenttime->tm_mon + 1) % 256,
            (currenttime->tm_mday) % 256,
            (currenttime->tm_hour) % 256,
            (currenttime->tm_min) % 256,
            (currenttime->tm_sec) % 256);
    fwrite(">>\nendobj\n", 10, 1, fout);
    fwrite("\n", 1, 1, fout);
    xref[xnum] = ftell(fout);
    xnum++;
    fwrite("xref\n0 ", 7, 1, fout);
    fprintf(fout, "%d\n%.10ld", xnum - 1, xref[0]);
    fwrite(" 65535 f\n", 9, 1, fout);
    for (kd = 1; kd < xnum - 1; kd++)
    {
        fprintf(fout, "%.10ld", xref[kd]);
        fwrite(" 00000 n\n", 9, 1, fout);
    }
    fwrite("trailer\n", 8, 1, fout);
    fwrite("\n<<\n/Info 11 0 R\n/Root 10 0 R\n/Size 12\n>>\nstartxref\n", 52, 1, fout);
    fprintf(fout, "%ld\n", xref[xnum - 1]);
    fwrite("%%EOF\n", 6, 1, fout);
    fclose(fout);

    if (buf)
        _TIFFfree(buf);
    if (bufbg)
        _TIFFfree(bufbg);
    if (buffg)
        _TIFFfree(buffg);
    return (0);
}
