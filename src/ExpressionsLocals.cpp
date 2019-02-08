#include <iostream>
#include <vec3.hpp>
#include <mat9.hpp>
#include <quat.hpp>

using namespace std;

int main (int argc, char const *argv[])
{
	double r = cos(M_PI/4.);
	vec3r n(0,r,r);
	vec3r t(0,-r,r);
	vec3r s(1,0,0);
	vec3r p0(0,0,0);
	mat9r P(n.x,n.y,n.z,t.x,t.y,t.z,s.x,s.y,s.z); // v_loc = P v_glob
	
	
	vec3r v_glob(0,-r,r);
	vec3r p_glob(0,0,0);
	quat Q_glob; //Q_glob.set_axis_angle(vec3r(1,0,0), M_PI/4.0);
	
	// calcul Q_loc
	quat qq; qq.set_rot_matrix(P.c_mtx());
	quat Q_loc = Q_glob.get_conjugated() * qq * Q_glob;
	
	// calcul v_loc
	vec3r v_loc = P * v_glob;
		
	// calcul p_loc
	vec3r p_loc = p_glob - p0;
	
	cout << "p_loc = " << p_loc << endl;
	cout << "v_loc = " << v_loc << endl;
	
	cout << "Q_glob * v_glob = " << Q_glob * v_glob << endl;
	cout << "Q_glob * v_loc = " << Q_glob * v_loc << endl;
	cout << "Q_loc * v_glob = " << Q_loc * v_glob << endl;
	
	//Q_glob * v_loc  =Q_glob * ( qq * v_glob)

	return 0;
}