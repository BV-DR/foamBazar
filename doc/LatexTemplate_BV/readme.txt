!---------------------------------!
!          Using bibTeX           !
!---------------------------------!


Instead of hard-coding the path to the .bib in .tex file, the following option can be given to bibtex :
 >> bibtex.exe --include-directory="C:\drsvn\drtools\LatexTemplate_BV" %

bibTeX can be edited manually but as it constantly growing, it is easier to use a bibTeX manager such as JabRef :
 --> https://www.jabref.org/ (works for Windows, Linux and MacOS)


!---------------------------------!
!     Template BV for LateX       !
!---------------------------------!

This is BV templates for Latex, based on the one made by Antoine with slight modifications to make it work with pdfLatex compilation.


1. Installation

	- On Windows:
	    -Copy "bv_article" et "bv_report" in the folder c:\..\MiKTeX 2.8\tex\latex
	    -Execute the "setting tools" of MiKTeX  (Menu Démarrer\Programmes\MiKTeX 2.8\Maintenance (Admin) )
	       -Refresh FNDB
	       -Update format
	    - proxy for the updates = euaproxy.bureauveritas.com:8080

	- On Linux:
	    - if not yet created, create folder ~/texmf/tex/latex
            - copy-paste "bv_article" and "bv_report" in ~/texmf/tex/latex


2.Usage

    -set the documentclass to bv_report or bv_article

    -Specific command of the BV templates :

       \titrerapport{Report}
       \nomaffaire{nom affaire}
       \noata{Numero ATA}
       \nont{Numero NT}
       \reva{0}\datea{21 Avril 2009}\writtenbya{author}\objreva{}\checkedbya{checker}
       \revb{1}\dateb{November 2009}\writtenbyb{author}\objrevb{}\checkedbyb{checker}
       \revc{1}\datec{November 2009}\writtenbyc{author}\objrevc{}\checkedbyc{checker}
       \marineconditions

       \begin{document}
       \bvpage

    -Works with classic LateX compilation as well as with pdfLateX
