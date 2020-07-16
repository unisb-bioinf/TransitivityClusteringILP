/** 
 * Copyright (C) 2020 Tim Kehl <tkehl@bioinf.uni-sb.de>
 *                    Kerstin Lenhof <klenhof@bioinf.uni-sb.de>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the Lesser GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * Lesser GNU General Public License for more details.
 *
 * You should have received a copy of the Lesser GNU General Public
 * License along with this program.
 * If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef TRANSITIVITY_CLUSTERING_ILP_CONFIG_MACROS_H
#define TRANSITIVITY_CLUSTERING_ILP_CONFIG_MACROS_H

#define GT2_EXPORT          __attribute__((visibility ("default")))
#define LOCAL           __attribute__((visibility ("hidden")))
#define EXTERN_VARIABLE extern __attribute__((visibility ("default")))

// This makro is needed to avoid 'unused parameter' warnings when using makros
#define _unused(x) ((void)(x))

#endif //TRANSITIVITY_CLUSTERING_ILP_CONFIG_MACROS_H
