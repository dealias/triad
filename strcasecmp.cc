#include <string.h>
/*
 * strcasecmp.c : Case-insensitive string comparison routines
 *
 * George Ferguson, ferguson@cs.rochester.edu, 23 Apr 1993.
 */

#define ISUPPER(C) ((C) >= 'A' && (C) <= 'Z')
#define TOLOWER(C) ((C) - 'A' + 'a')
#define NORMAL(C) (ISUPPER(C) ? TOLOWER(C) : (C))

int strcasecmp(const char *s1, const char *s2)
{
    char c1,c2;

    while (1) {
		c1 = NORMAL(*s1);
		c2 = NORMAL(*s2);
		if (c1 != c2) return(c1-c2);
		else if (c1 == '\0') return(0);
		else {
			s1++;
			s2++;
		}
    }
    /* NOT REACHED */
}

int strncasecmp(const char *s1, const char *s2, size_t n)
{
    char c1,c2;

    while (n-- > 0) {
		c1 = NORMAL(*s1);
		c2 = NORMAL(*s2);
		if (c1 != c2) return(c1-c2);
		else if (c1 == '\0') return(0);
		else {
			s1++;
			s2++;
		}
    }
    return(0);
}
