/*******************************************************************************
 *  @file edgebounds_obj.h
 *  @brief EDGEBOUNDS Object functions.
 *
 *  @author Dave Rich (devrek)
 *  @bug Lots.
 *******************************************************************************/


#ifndef _EDGEBOUNDS_OBJ_H
#define _EDGEBOUNDS_OBJ_H

void edgebounds_Create(EDGEBOUNDS *edg);

void edgebounds_Destroy(EDGEBOUNDS *edg);

void edgebounds_Resize(EDGEBOUNDS *edg);

void edgebounds_Print(EDGEBOUNDS *edg);

void edgebounds_Merge(EDGEBOUNDS *edg_fwd,
                     EDGEBOUNDS *edg_bck,
                     EDGEBOUNDS *edg_new);

void edgebounds_Reorient(EDGEBOUNDS *edg);

#endif /* _EDGEBOUNDS_OBJ_H */