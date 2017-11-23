#include "quad.h"

// Create a square quadrant with a given length and midpoints (xmid,ymid)
Quad::Quad(double xmid, double ymid, double length) {
	this->xmid = xmid;
	this->ymid = ymid;
	this->length = length;
}

// How long is this quadrant?
double Quad::get_length() {
	return this->length;
}

// Check if the current quadrant contains a point
bool Quad::contains(double xmid, double ymid) {
	if (xmid <= this->xmid + this->length / 2.0 && xmid >= this->xmid - this->length / 2.0 && ymid <= this->ymid + this->length / 2.0 && ymid >= this->ymid - this->length / 2.0) {
		return true;
	}
	else {
		return false;
	}
}

// Create subdivisions of the current quadrant
// Subdivision: Northwest quadrant
Quad Quad::NW() {
	Quad newquad = Quad(this->xmid - this->length / 4.0, this->ymid + this->length / 4.0, this->length / 2.0);
	return newquad;
}

// Subdivision: Northeast quadrant
Quad Quad::NE() {
	Quad newquad = Quad(this->xmid + this->length / 4.0, this->ymid + this->length / 4.0, this->length / 2.0);
	return newquad;
}

// Subdivision: Southwest quadrant
Quad Quad::SW() {
	Quad newquad = Quad(this->xmid - this->length / 4.0, this->ymid - this->length / 4.0, this->length / 2.0);
	return newquad;
}

// Subdivision: Southeast quadrant
Quad Quad::SE() {
	Quad newquad = Quad(this->xmid + this->length / 4.0, this->ymid - this->length / 4.0, this->length / 2.0);
	return newquad;
}