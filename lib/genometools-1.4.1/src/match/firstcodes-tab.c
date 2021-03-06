/*
  Copyright (c) 2011 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#ifndef S_SPLINT_S
#include "core/xansi_api.h"
#endif
#include "core/fa.h"
#include "core/assert_api.h"
#include "core/divmodmul.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#ifdef SKDEBIG
#include "core/disc_distri.h"
#endif
#include "core/log.h"
#include "core/spacecalc.h"
#include "core/hashmap-generic.h"
#include "core/arraydef.h"
#include "firstcodes-tab.h"
#include "firstcodes-psbuf.h"

static void gt_firstcodes_countocc_new(GtFirstcodesspacelog *fcsl,
                                       GtFirstcodestab *fct,
                                       unsigned long numofsequences)
{
  fct->countocc_small = gt_malloc((size_t) (numofsequences+1) *
                                  sizeof (*fct->countocc_small));
  GT_FCI_ADDWORKSPACE(fcsl,"countocc_small",
                      sizeof (*fct->countocc_small) *
                      (numofsequences+1));
  fct->countocc_exceptions = ul_u32_gt_hashmap_new();
  gt_assert(fct->countocc_exceptions != NULL);
  fct->outfilenameleftborder = NULL;
  fct->leftborder_samples = NULL;
#ifdef _LP64
  fct->modvaluebits = 32U; /* XXX remove the following later */
  if (fct->modvaluebits == 32U)
  {
    fct->modvaluemask = UINT32_MAX;
  } else
  {
    fct->modvaluemask = (uint32_t) ((1UL << fct->modvaluebits) - 1);
  }
  GT_INITARRAY(&fct->bitchangepoints,GtUlong);
#endif
}

static void gt_firstcodes_countocc_resize(GtFirstcodesspacelog *fcsl,
                                          GtFirstcodestab *fct,
                                          unsigned long numofdifferentcodes)
{
  fct->countocc_small = gt_realloc(fct->countocc_small,
                                   sizeof (*fct->countocc_small) *
                                           (numofdifferentcodes+1));
  GT_FCI_SUBTRACTADDWORKSPACE(fcsl,"countocc_small",
                              sizeof (*fct->countocc_small) *
                              (numofdifferentcodes+1));
}

#ifdef SKDEBUG

typedef struct
{
  unsigned long smallcount, smallsum,
                largecount, largesum,
                hugecount, hugesum;
} GtCountdistri_info;

static void gt_firstcodes_evaluate_distvalue(unsigned long key,
                                             unsigned long long value,
                                             void *data)
{
  GtCountdistri_info *cdi = (GtCountdistri_info *) data;

  gt_assert(value > 0);
  if (key <= UINT8_MAX)
  {
    cdi->smallcount++;
    cdi->smallsum += value;
  } else
  {
    if (key <= UINT16_MAX)
    {
      cdi->largecount++;
      cdi->largesum += value;
    } else
    {
      cdi->hugecount++;
      cdi->hugesum += value;
    }
  }
}

static void gt_firstcodes_evaluate_countdistri(const GtDiscDistri *countdistri)
{
  GtCountdistri_info cdi;
  unsigned long sum;
  size_t spacenow, spacedirectstore, spaceopt, spacewithhash;;

  cdi.smallcount = 0;
  cdi.smallsum = 0;
  cdi.largecount = 0;
  cdi.largesum = 0;
  cdi.hugecount = 0;
  cdi.hugesum = 0;
  gt_disc_distri_foreach(countdistri,gt_firstcodes_evaluate_distvalue,&cdi);
  sum = cdi.smallsum + cdi.largesum + cdi.hugesum;
  gt_log_log("small=%lu,%lu (%.2f)",cdi.smallcount,cdi.smallsum,
          (double) cdi.smallsum/sum);
  gt_log_log("large=%lu,%lu (%.2f)",cdi.largecount,cdi.largesum,
          (double) cdi.largesum/sum);
  gt_log_log("huge=%lu,%lu (%.2f)",cdi.hugecount,cdi.hugesum,
          (double) cdi.largesum/sum);
  spacenow = sizeof (uint32_t) * sum;
  spaceopt = sizeof (uint8_t) * sum;
  spacedirectstore = sizeof (uint32_t) * cdi.largesum;
  spacewithhash = (2 * sizeof (void *)) * cdi.largecount;
  gt_log_log("spacenow=%.2f, spaceopt (direct)=%2.f spaceopt (hash) =%.2f",
          GT_MEGABYTES(spacenow),
          GT_MEGABYTES(spaceopt+spacedirectstore),
          GT_MEGABYTES(spaceopt+spacewithhash));
}
#endif

#ifdef SKDEBUG
static void checkcodesorder(const unsigned long *tab,unsigned long len,
                            bool allowequal)
{
  unsigned long idx;

  for (idx=1UL; idx < len; idx++)
  {
    gt_assert(tab[idx-1] < tab[idx] || (allowequal && tab[idx-1] == tab[idx]));
  }
}
#endif

unsigned long gt_firstcodes_remdups(unsigned long *allfirstcodes,
                                    GtFirstcodesspacelog *fcsl,
                                    GtFirstcodestab *fct,
                                    unsigned long numofsequences,
                                    Gtmarksubstring *markprefix,
                                    Gtmarksubstring *marksuffix,
                                    GtLogger *logger)
{
  if (numofsequences == 0)
  {
    fct->differentcodes = 0;
  } else
  {
    unsigned long numofdifferentcodes, *storeptr, *readptr;
    bool firstincrement;

    gt_firstcodes_countocc_new(fcsl,fct,numofsequences);
    gt_firstcodes_countocc_increment(fct,0,true); /* first increment */
    gt_marksubstring_mark(markprefix,allfirstcodes[0]);
    gt_marksubstring_mark(marksuffix,allfirstcodes[0]);
    for (storeptr = allfirstcodes, readptr = allfirstcodes+1;
         readptr < allfirstcodes + numofsequences;
         readptr++)
    {
      if (*storeptr != *readptr)
      {
        storeptr++;
        *storeptr = *readptr;
        firstincrement = true;
      } else
      {
        firstincrement = false;
      }
      gt_firstcodes_countocc_increment(fct,(unsigned long)
                                       (storeptr - allfirstcodes),
                                       firstincrement);
      gt_marksubstring_mark(markprefix,*readptr);
      gt_marksubstring_mark(marksuffix,*readptr);
    }
    numofdifferentcodes = (unsigned long) (storeptr - allfirstcodes + 1);
    if (numofdifferentcodes < numofsequences)
    {
      /* reduce the memory requirement, as the duplicated elements are not
         needed */
      gt_firstcodes_countocc_resize(fcsl,fct,numofdifferentcodes);
#ifdef SKDEBUG
      checkcodesorder(allfirstcodes,numofdifferentcodes,false);
#endif
    }
    fct->differentcodes = numofdifferentcodes;
  }
  gt_logger_log(logger,"number of different first codes=%lu (%.2f%%) "
                       "in %lu sequences",
                fct->differentcodes,
                100.00 * (double) fct->differentcodes/numofsequences,
                numofsequences);
  return fct->differentcodes;
}

static uint32_t gt_firstcodes_countocc_get(const GtFirstcodestab *fct,
                                           unsigned long idx)
{
  if (fct->countocc_small[idx] > 0)
  {
    return (uint32_t) fct->countocc_small[idx];
  } else
  {
    uint32_t *valueptr = ul_u32_gt_hashmap_get(fct->countocc_exceptions,idx);

    gt_assert(valueptr != NULL);
    return *valueptr + (uint32_t) GT_FIRSTCODES_MAXSMALL;
  }
}

#define GT_PARTIALSUM_COUNT_GET(IDX)   gt_firstcodes_countocc_get(fct,IDX)

#ifdef _LP64
#define GT_PARTIALSUM_LEFTBORDER_SET(BUF,VALUE)\
        if ((BUF)->nextfree == (BUF)->allocated)\
        {\
          gt_leftborderbuffer_flush(BUF);\
        }\
        (BUF)->spaceuint32_t[(BUF)->nextfree++]\
          = (uint32_t) ((VALUE) & fct->modvaluemask)
#else
#define GT_PARTIALSUM_LEFTBORDER_SET(BUF,VALUE)\
        if ((BUF)->nextfree == (BUF)->allocated)\
        {\
          gt_leftborderbuffer_flush(BUF);\
        }\
        (BUF)->spaceuint32_t[(BUF)->nextfree++] = (uint32_t) (VALUE)
#endif

#define GT_FIRSTCODES_ADD_SAMPLE(PARTSUM)\
        gt_assert(samplecount < fct->numofsamples);\
        fct->leftborder_samples[samplecount++] = PARTSUM

unsigned long gt_firstcodes_partialsums(GtFirstcodesspacelog *fcsl,
                                        GtFirstcodestab *fct,
                                        GT_UNUSED unsigned long
                                                            expectedlastpartsum)
{
  unsigned long idx, partsum, maxbucketsize, bitmask, samplecount = 0,
                spacewithhashmap = 0, spacewithouthashmap = 0;
  uint32_t currentcount;
  GtLeftborderOutbuffer *leftborderbuffer_all = NULL;
#ifdef _LP64
  const unsigned int btp = gt_determinebitspervalue(expectedlastpartsum);
  unsigned long exceedvalue = 1UL << fct->modvaluebits;
#endif
#ifdef SKDEBUG
  GtDiscDistri *countdistri = gt_disc_distri_new();
#endif

  gt_assert(fct->differentcodes < UINT32_MAX);
  gt_log_log("hashmap_addcount=%lu (%.2f%%)",fct->hashmap_addcount,
                  100.0 * (double) fct->hashmap_addcount/
                                   fct->differentcodes);
  gt_log_log("hashmap_incrementcount=%lu (%.2f%%)",
                  fct->hashmap_incrementcount,
                  100.0 * (double) fct->hashmap_incrementcount/
                                   fct->all_incrementcount);
  gt_log_log("hashmap_getcount=%lu (%.2f%%)",
                  fct->hashmap_getcount,
                  100.0 * (double) fct->hashmap_getcount/
                                   fct->all_incrementcount);

#ifdef _LP64
  if (btp <= fct->modvaluebits)
  {
    fct->bitchangepoints.allocatedGtUlong = 0;
    fct->bitchangepoints.spaceGtUlong = NULL;
  } else
  {
    fct->bitchangepoints.allocatedGtUlong = 1UL << (btp - fct->modvaluebits);
    gt_log_log("lastpartsum=%lu, bitchangepoints.allocated=%lu",
              expectedlastpartsum,fct->bitchangepoints.allocatedGtUlong);
    fct->bitchangepoints.spaceGtUlong
      = gt_malloc(sizeof (*fct->bitchangepoints.spaceGtUlong)
                  * fct->bitchangepoints.allocatedGtUlong);
  }
  fct->bitchangepoints.nextfreeGtUlong = 0;
#endif
  currentcount = GT_PARTIALSUM_COUNT_GET(0);
  partsum = (unsigned long) currentcount;
  maxbucketsize = (unsigned long) currentcount;
#ifdef SKDEBUG
  gt_disc_distri_add(countdistri,(unsigned long) currentcount);
#endif
  fct->sampleshift = 9U;
  while (true)
  {
    fct->sampledistance = 1UL << fct->sampleshift;
    if (fct->sampledistance < fct->differentcodes)
    {
      break;
    }
    fct->sampleshift--;
  }
  bitmask = fct->sampledistance - 1;
  fct->numofsamples = 1UL + 1UL + fct->differentcodes/fct->sampledistance;
  fct->leftborder_samples = gt_malloc(sizeof (*fct->leftborder_samples) *
                                      fct->numofsamples);
  GT_FCI_ADDWORKSPACE(fcsl,"leftborder_samples",
                      sizeof (*fct->leftborder_samples) * fct->numofsamples);
  GT_FIRSTCODES_ADD_SAMPLE(partsum);
  leftborderbuffer_all = gt_leftborderbuffer_new("leftborder",fcsl);
  GT_PARTIALSUM_LEFTBORDER_SET(leftborderbuffer_all,partsum);
  for (idx = 1UL; idx < fct->differentcodes; idx++)
  {
    currentcount = GT_PARTIALSUM_COUNT_GET(idx);
#ifdef _LP64
    gt_assert(currentcount <= fct->modvaluemask);
#endif
#ifdef SKDEBUG
    gt_disc_distri_add(countdistri,(unsigned long) currentcount);
#endif
    if (maxbucketsize < (unsigned long) currentcount)
    {
      maxbucketsize = (unsigned long) currentcount;
    }
    partsum += currentcount;
#ifdef _LP64
    if (fct->bitchangepoints.allocatedGtUlong > 0 && partsum >= exceedvalue)
    {
      gt_assert(idx > 0 && fct->bitchangepoints.nextfreeGtUlong <
                           fct->bitchangepoints.allocatedGtUlong);
      gt_assert(fct->bitchangepoints.spaceGtUlong != NULL);
      fct->bitchangepoints.spaceGtUlong
        [fct->bitchangepoints.nextfreeGtUlong++] = idx-1;
      exceedvalue
        = ((exceedvalue >> fct->modvaluebits) + 1) << fct->modvaluebits;
    }
#endif
    if ((idx & bitmask) == 0)
    {
      GT_FIRSTCODES_ADD_SAMPLE(partsum);
    }
    GT_PARTIALSUM_LEFTBORDER_SET(leftborderbuffer_all,partsum);
  }
  GT_PARTIALSUM_LEFTBORDER_SET(leftborderbuffer_all,partsum);
  fct->outfilenameleftborder
      = gt_leftborderbuffer_delete(leftborderbuffer_all,fcsl,
                                   gt_firstcodes_leftborder_entries(fct));
  if (partsum > fct->leftborder_samples[samplecount-1])
  {
    GT_FIRSTCODES_ADD_SAMPLE(partsum);
  } else
  {
    gt_assert(partsum == fct->leftborder_samples[samplecount-1]);
  }
  gt_assert(expectedlastpartsum == partsum);
  fct->numofsamples = samplecount-1;
#ifdef SKDEBUG
  gt_firstcodes_evaluate_countdistri(countdistri);
  gt_disc_distri_delete(countdistri);
#endif
  gt_assert (fct->countocc_small != NULL);
  gt_free(fct->countocc_small);
  GT_FCI_SUBTRACTWORKSPACE(fcsl,"countocc_small");
  fct->countocc_small = NULL;
  if (fct->hashmap_addcount > 0 && gt_ma_bookkeeping_enabled())
  {
    spacewithhashmap = gt_ma_get_space_current() + gt_fa_get_space_current();
  }
  gt_hashtable_delete(fct->countocc_exceptions);
  if (fct->hashmap_addcount > 0 && gt_ma_bookkeeping_enabled())
  {
    unsigned long hashmapspace;

    spacewithouthashmap = gt_ma_get_space_current() + gt_fa_get_space_current();
    gt_assert(spacewithouthashmap < spacewithhashmap);
    hashmapspace = spacewithhashmap - spacewithouthashmap;
    gt_log_log("space for hashmap=%.2f (%lu bytes per entry)",
               GT_MEGABYTES(hashmapspace),hashmapspace/fct->hashmap_addcount);
  }
  fct->countocc_exceptions = NULL;
  return maxbucketsize;
}

unsigned long gt_firstcodes_get_sample(const GtFirstcodestab *fct,
                                       unsigned long idx)
{
  gt_assert(idx <= fct->numofsamples);
  return fct->leftborder_samples[idx];
}

unsigned long gt_firstcodes_get_leftborder(const GtFirstcodestab *fct,
                                           unsigned long idx)
{
#ifdef _LP64
  GT_CHANGEPOINT_GET(changepoint);

  return (unsigned long) fct->leftborder[idx]
                         + (changepoint << fct->modvaluebits);
#else
  return (unsigned long) fct->leftborder[idx];
#endif
}

unsigned long gt_firstcodes_leftborder_entries(const GtFirstcodestab *fct)
{
  return fct->differentcodes + 1;
}

unsigned long gt_firstcodes_numofsamples(const GtFirstcodestab *fct)
{
  return fct->numofsamples;
}

unsigned long gt_firstcodes_findfirstsamplelarger(const GtFirstcodestab *fct,
                                                  unsigned long suftaboffset)
{
  unsigned long left = 0, right, mid, midval, found;

  right = found = fct->numofsamples;
  while (left+1 < right)
  {
    mid = GT_DIV2(left+right);
    midval = gt_firstcodes_get_sample(fct,mid);
    if (suftaboffset == midval)
    {
      return mid;
    }
    if (suftaboffset < midval)
    {
      found = mid;
      right = mid - 1;
    } else
    {
      left = mid + 1;
    }
  }
  gt_assert(suftaboffset <= gt_firstcodes_get_sample(fct,found));
  return found;
}

unsigned long gt_firstcodes_sample2full(const GtFirstcodestab *fct,
                                        unsigned long idx)
{
  gt_assert(idx <= fct->numofsamples);
  if (idx < fct->numofsamples)
  {
    return idx << fct->sampleshift;
  }
  return fct->differentcodes - 1;
}

void gt_firstcodes_samples_delete(GtFirstcodesspacelog *fcsl,
                                  GtFirstcodestab *fct)
{
  if (fct->leftborder_samples != NULL)
  {
    gt_free(fct->leftborder_samples);
    GT_FCI_SUBTRACTWORKSPACE(fcsl,"leftborder_samples");
    fct->leftborder_samples = NULL;
  }
}

void gt_firstcodes_countocc_delete(GtFirstcodesspacelog *fcsl,
                                   GtFirstcodestab *fct)
{
  if (fct->countocc_small != NULL)
  {
    GT_FCI_SUBTRACTWORKSPACE(fcsl,"countocc_small");
    gt_free(fct->countocc_small);
    fct->countocc_small = NULL;
  }
  gt_hashtable_delete(fct->countocc_exceptions);
  fct->countocc_exceptions = NULL;
}

void gt_firstcodes_tab_delete(GtFirstcodesspacelog *fcsl,GtFirstcodestab *fct)
{
  gt_firstcodes_samples_delete(fcsl,fct);
  gt_str_delete(fct->outfilenameleftborder);
  fct->outfilenameleftborder = NULL;
#ifdef _LP64
  GT_FREEARRAY(&fct->bitchangepoints,GtUlong);
#endif
}

void gt_firstcodes_countocc_setnull(GtFirstcodestab *fct)
{
  fct->leftborder = NULL;
  fct->countocc_small = NULL;
  fct->leftborder_samples = NULL;
  fct->countocc_exceptions = NULL;
  fct->differentcodes = 0;
  fct->lastincremented_idx = 0;
  fct->lastincremented_valueptr = NULL;
  fct->hashmap_addcount = 0;
  fct->hashmap_incrementcount = 0;
  fct->all_incrementcount = 0;
  fct->hashmap_getcount = 0;
  fct->outfilenameleftborder = NULL;
#ifdef _LP64
  GT_INITARRAY(&fct->bitchangepoints,GtUlong);
#endif
}

uint32_t **gt_firstcodes_leftborder_address(GtFirstcodestab *fct)
{
  return &fct->leftborder;
}

void gt_firstcodes_leftborder_remap(GtFirstcodestab *fct,uint32_t *ptr)
{
  fct->leftborder = ptr;
}

const GtStr *gt_firstcodes_outfilenameleftborder(const GtFirstcodestab *fct)
{
  return fct->outfilenameleftborder;
}
