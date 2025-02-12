/*
 * VK-Gong code 
 *
 * Author of this version: Àngels Aragonès
 *
 * Based on VK-Gong for Matlab and NeSS Framework Code (Copyright (c) The University of Edinburgh, 2014-2015. All rights reserved. Author: James Perry)
 *
 * The logging system.
 */
#ifndef _LOGGER_H_
#define _LOGGER_H_

// logs a message (using printf-style arguments), this may go to stderr or to a file
// depending on environment settings
void logMessage(int level, char *msg, ...);

void logClose();

#endif
