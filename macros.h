/**
 * Transitivity Clustering ILP
 * Copyright (C) 2019 Tim Kehl <tkehl@bioinf.uni-sb.de>
 */

#ifndef TRANSITIVITY_CLUSTERING_ILP_CONFIG_MACROS_H
#define TRANSITIVITY_CLUSTERING_ILP_CONFIG_MACROS_H

#define GT2_EXPORT          __attribute__((visibility ("default")))
#define LOCAL           __attribute__((visibility ("hidden")))
#define EXTERN_VARIABLE extern __attribute__((visibility ("default")))

// This makro is needed to avoid 'unused parameter' warnings when using makros
#define _unused(x) ((void)(x))

#endif //TRANSITIVITY_CLUSTERING_ILP_CONFIG_MACROS_H
