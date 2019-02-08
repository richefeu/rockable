#ifndef PSIO_HPP_72530F82
#define PSIO_HPP_72530F82

struct psio {
	FILE *psout;
	float sca;

	psio(const char * nom, float xb1, float yb1, float xb2, float yb2 , float SCA = 1.0, float vesca = 5.0): sca(SCA) {
		psout = fopen (nom, "w");
		if (psout == NULL) {
			fprintf(stderr, "Impossible de creer le fichier %s\n", nom);
		}

		fprintf(psout, "%%!PS-Adobe-3.0 EPSF-3.0\n");
		fprintf(psout, "%%%%BoundingBox: %g %g %g %g\n", floor(sca * xb1), floor(sca * yb1), ceil(sca * xb2), ceil(sca * yb2));
		fprintf(psout, "%%%%HiResBoundingBox: %g %g %g %g\n", sca * xb1, sca * yb1, sca * xb2, sca * yb2);
		fprintf(psout, "%%%%EndComments\n");

		fprintf(psout, "1 setlinecap\n");

		// Procedure cprint: texte centre >> (texte) x y cprint
		fprintf(psout, "/cprint{moveto dup stringwidth pop 2 div neg 0 rmoveto show}def\n");
		fprintf(psout, "/rprint{moveto dup stringwidth pop 1 div neg 0 rmoveto show}def\n");
		fprintf(psout, "/Times-Roman findfont 15 scalefont setfont\n");

		// Procedure C: trace un cercle de centre xy de rayon r et niveau de gris g (1=blanc, 0=noir) >> x y r g C
		fprintf(psout, "/C{newpath 4 1 roll 0 360 arc gsave setgray fill grestore stroke}def\n");

		// Procedure Arc: trace un arc de cercle de centre xy de rayon r et niveau de gris g (1=blanc, 0=noir)
		// de l'angle deb a fin >> x y r deb fin g Arc
		fprintf(psout, "/Arc{newpath 6 1 roll arc gsave setgray fill grestore stroke}def\n");

		// Procedure Lic: ligne intercentre de largeur donnee largeur t,
		// du point 0 au point 1 >> x0 y0 x1 y1 t Lic
		//fprintf(psout,"/Lic{gsave setrgbcolor 0.1 mul setlinewidth newpath moveto lineto stroke grestore}def\n");
		//fprintf(psout,"/Lic{gsave 0.1 mul setlinewidth newpath moveto lineto stroke grestore}def\n");
		fprintf(psout, "/Lic{gsave setlinewidth newpath moveto lineto stroke grestore}def\n");

		// Procedure Lis: ligne intercentre de largeur donnee largeur fixe,
		// du point 0 au point 1 >> x0 y0 x1 y1 Lis
		fprintf(psout, "/Lis{gsave newpath moveto lineto stroke grestore}def");

		// Procedure Ve: Vecteur >> longueur angle x0 y0 Ve
		fprintf(psout, "/fref{2 setlinewidth newpath 0 0 moveto 86 0 lineto stroke\n");
		fprintf(psout, " newpath 86 0 moveto 80 -6 lineto 100 0 lineto 80 6 lineto closepath fill}def\n");
		fprintf(psout, "/Ve{gsave translate rotate vesca mul dup scale fref grestore}def\n");
		fprintf(psout, "/vesca %g def\n", vesca);

		// util
		fprintf(psout, "/m {moveto} def\n");
		fprintf(psout, "/l {lineto} def\n");
		fprintf(psout, "/s {stroke} def\n");
		fprintf(psout, "/n {newpath} def\n");
		fprintf(psout, "/c {closepath} def\n");
	}

	~psio() {
		fclose(psout);
	}

	void setrgbcolor(float r, float g, float b) {
		fprintf(psout, "%g %g %g setrgbcolor\n", r, g, b);
	}

	void raw(const char * txt) {
		fprintf(psout, "%s\n", txt);
	}

	void line(float x1, float y1, float x2, float y2) {
		fprintf(psout, "%g %g %g %g Lis\n", sca * x1, sca * y1, sca * x2, sca * y2);
	}

	void inter_center(float x1, float y1, float x2, float y2, float t) {
		fprintf(psout, "%g %g %g %g %g Lic\n", sca * x1, sca * y1, sca * x2, sca * y2, sca * t);
	}

	void vector(float l, float ang, float x, float y) {
		fprintf(psout, "%g %g %g %g Ve\n", sca * x, sca * y, ang, sca * l);
	}

	void circle(float x, float y, float R, float gray) {
		fprintf(psout, "%g %g %g %g C\n", sca * x, sca * y, sca * R, gray);
	}

	void arc(float x, float y, float R, float deb, float fin) {
		fprintf(psout, "%g %g %g %g %g 1 Arc\n", sca * x, sca * y, sca * R, deb, fin);
	}
};


#endif /* end of include guard: PSIO_HPP_72530F82 */

