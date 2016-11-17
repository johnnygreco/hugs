from .. import SynthFactory
from ...utils import get_test_exp, remove_mask_planes

def test_SynthFactory():
    size = 500
    sf = SynthFactory(num_synths=10,
                      pset_lims={'X0': [0, size], 'Y0': [0, size]})
    ps = sf.get_psets()
    assert len(ps)==10

    img = sf.create_image((size, size))
    assert img.shape==(size, size)

    exp = get_test_exp()
    msk = exp.getMaskedImage().getMask()
    remove_mask_planes(msk, ['CR'])
    sf.inject(exp)
    assert 'SYNTH' in msk.getMaskPlaneDict().keys()
    remove_mask_planes(msk, ['SYNTH'])
