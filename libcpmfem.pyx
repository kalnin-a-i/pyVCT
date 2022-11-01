cdef extern from "libcpmfem.h":
    int cpmfem(int NCX, int NCY, 
	double PART,
	double GN_CM,
	double GN_FB,
	double TARGETVOLUME_CM,
	double DETACH_CM,
	double DETACH_FB,
	double INELASTICITY_FB,
	double INELASTICITY_CM,
	double JCMMD,
	double JFBMD,
	double JCMCM,
	double JFBFB,
	double JFBCM,
	double UNLEASH_CM,
	double UNLEASH_FB,
	double LMAX_CM,
	double LMAX_FB,
	double MAX_FOCALS_CM,
	double MAX_FOCALS_FB)

def py_cpmfem(NCX, NCY, PART, GN_CM, GN_FB, TARGETVOLUME_CM, DETACH_CM, DETACH_FB, INELASTICITY_CM,
              INELASTICITY_FB, JCMMD, JFBMD, JCMCM, JFBFB, JFBCM, UNLEASH_CM, UNLEASH_FB, LMAX_CM, LMAX_FB,
              MAX_FOCALS_CM, MAX_FOCALS_FB):
    return cpmfem(NCX, NCY, PART, GN_CM, GN_FB, TARGETVOLUME_CM, DETACH_CM, DETACH_FB, INELASTICITY_CM,
              INELASTICITY_FB, JCMMD, JFBMD, JCMCM, JFBFB, JFBCM, UNLEASH_CM, UNLEASH_FB, LMAX_CM, LMAX_FB,
              MAX_FOCALS_CM, MAX_FOCALS_FB)