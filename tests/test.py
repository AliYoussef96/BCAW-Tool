import pytest
from BCAWT import BCAWT
from BCAWT import ATCG3
from BCAWT import CA_RSCU
from BCAWT import GC123
from BCAWT import GRAVY_AROMO
from BCAWT import PR2_plot_data
from BCAWT import P2_index


def test_bad_args_BCAWT():
    # make sure bad arguments raise errors
    with pytest.raises(TypeError):
        BCAWT.BCAW('Ecoli','test_demo',fasta = True,txt = True, Auto=True) #fasta and txt = True
    with pytest.raises(TypeError):
        BCAWT.BCAW('Ecoli','test_demo',fasta = False,txt = False, Auto=True) #fasta and txt = False
    with pytest.raises(TypeError):
        BCAWT.BCAW('Ecoli','test_demo',Auto=True) #fasta and txt = False by default
    with pytest.raises(TypeError):
        BCAWT.BCAW('Ecoli','test_demo',"Ecoli",fasta = True,txt = True, Auto=True) #Both genes reference set and Auto are specified

#make sure other modules return the true result or not

def test_ATCG3_result():
    assert ATCG3.ACTG3("NNANNTNNGNNC", A=True) == 25.0
    assert ATCG3.ACTG3("NNANNTNNGNNC", T=True) == 25.0
    assert ATCG3.ACTG3("NNANNTNNGNNC", C=True) == 25.0
    assert ATCG3.ACTG3("NNANNTNNGNNC", G=True) == 25.0

def test_cu_rscu_result():
    df = CA_RSCU.CA_RSCU("GGGCCCAAATTT","test_gene_name")
    assert list(df.index) == ['AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATT', 'CAA', 'CAC', 'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 'GTC', 'GTG', 'GTT', 'TAC', 'TAT', 'TCA', 'TCC', 'TCG', 'TCT', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 'TTT']

    assert list(df.test_gene_name) == [2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0]

def test_gc123_result():
    assert GC123.GC1("GGGCCCAAATTT") == 50.0
    assert GC123.GC2("GGGCCCAAATTT") == 50.0
    assert GC123.GC3("GGGCCCAAATTT") == 50.0
    assert GC123.GC12("GGGCCCAAATTT") == 50

def test_GRAVY_AROMO_result():
    assert GRAVY_AROMO.GRAvy_ARomo("GGGCCCAAATTT",G=True) == -0.7750000000000001
    assert GRAVY_AROMO.GRAvy_ARomo("GGGCCCAAATTT",A=True) == 0.25

def test_P2_index_result():
    assert P2_index.P2_index("AACGGC" ,p2 = True ) == 0.5
    assert P2_index.P2_index("AACGGC" ,wwc = True ) == 1
    assert P2_index.P2_index("AACGGC" ,sst = True ) == 0
    assert P2_index.P2_index("AACGGC" ,wwy = True ) == 1
    assert P2_index.P2_index("AACGGC" ,ssy = True ) == 1

def test_pr2_plot_result():
    assert PR2_plot_data.PR2_plot("NNANNTNNGNNC" ,o = True ) == 0.5
    assert PR2_plot_data.PR2_plot("NNANNTNNGNNC" ,a = True ) == 0.5
