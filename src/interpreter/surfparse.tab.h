#ifndef BISON_SURFPARSE_TAB_H
# define BISON_SURFPARSE_TAB_H

# ifndef YYSTYPE
#  define YYSTYPE int
#  define YYSTYPE_IS_TRIVIAL 1
# endif
# define	STANDARD_COMMAND	257
# define	IDENTIFIER	258
# define	STRING	259
# define	INTEGER	260
# define	REAL	261


extern YYSTYPE yylval;

#endif /* not BISON_SURFPARSE_TAB_H */
