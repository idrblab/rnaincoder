from padelpy import padeldescriptor
import pandas as pd
import staticsMethods as statics
import Dict_methods as Dict_methods
dictMe = {
    "Autocorrelation descriptors (1D)": "a1",
    "Barysz matrix descriptors (1D)": "a2",
    "Constitutional descriptors (1D)": "a3",
    "Topology related descriptors (1D)": "a4",
    "Physicochemical descriptors (1D)": "a5",
    "3D autocorrelation based (1D)": "b1",
    "Charged partial surface area (1D)": "b2",
    "Gravitational index based (1D)": "b3",
    "Length over breadth based (1D)": "b4",
    "Moment inertia descriptors (1D)": "b5",
    "Petitjean shape index based (1D)": "b6",
    "Radial distribution function (1D)": "b7",
    "Pubchem fingerprint (1D)": "c1"
}
def Moleculer_methods(uploadsdf_path,ResfilePath,InterFilepath,methodslist):
    padeldescriptor(mol_dir=uploadsdf_path,
                            d_file=ResfilePath,
                            convert3d=False,
                            retain3d=False,
                            retainorder=True,
                            d_2d=True,
                            d_3d=True,
                            fingerprints=True,
                            sp_timeout=60)

    indexa1 = []
    indexa2 = []
    indexa3 = []
    indexa4 = []
    indexa5 = []
    indexb1 = []
    indexb2 = []
    indexb3 = []
    indexb4 = []
    indexb5 = []
    indexb6 = []
    indexb7 = []
    indexc1 = []
    for method in methodslist:
        if dictMe[method] == 'a1':
            indexa1 = list(range(22,368,1))
        if dictMe[method] == 'a2':
            indexa2 = list(range(368,459,1))
        if dictMe[method] == 'a3':
            indexa3 = list(range(638,650,1))
        if dictMe[method] == 'a4':
            indexa4 = list(range(6,22,1))+list(range(459,476,1))+list(range(477,638,1))+list(range(652,1253,1))+list(range(1255,1274,1))+list(range(1280,1375,1))+list(range(1376,1412,1))+list(range(1414,1434,1))
        if dictMe[method] == 'a5':
            indexa5 = list(range(1,6,1))+list(range(476,477,1))+list(range(650,652,1))+list(range(1253,1255,1))+list(range(1274,1280,1))+list(range(1375,1376,1))+list(range(1412,1414,1))+list(range(1434,1445,1))
        if dictMe[method] == 'b1':
            indexb1 = [1445,1446,1447,1448,1449,1450,1451,1452,1453,1454,1455,1456,1457,1458,1459,1460,1461,1462,1463,1464,1465,1466,1467,1468,1469,1470,1471,1472,1473,1474,1475,1476,1477,1478,1479,1480,1481,1482,1483,1484,1485,1486,1487,1488,1489,1490,1491,1492,1493,1494,1495,1496,1497,1498,1499,1500,1501,1502,1503,1504,1505,1506,1507,1508,1509,1510,1511,1512,1513,1514,1515,1516,1517,1518,1519,1520,1521,1522,1523,1524,1785,1786,1787,1788,1789,1790,1791,1792,1793,1794,1795,1796,1797,1798,1799,1800,1801,1802,1803,1804,1805,1806,1807,1808,1809,1810,1811,1812,1813,1814,1815,1816,1817,1818,1819,1820,1821,1822,1823,1824,1825,1826,1827,1828,1829,1830,1831,1832,1833,1834,1835,1836,1837,1838,1839,1840,1841,1842,1843,1844,1845,1846,1847,1848,1849,1850,1851,1852,1853,1854,1855,1856,1857,1858,1859,1860,1861,1862,1863,1864,1865,1866,1867,1868,1869,1870,1871,1872,1873,1874,1875]
        if dictMe[method] == 'b2':
            indexb2 = [1525,1526,1527,1528,1529,1530,1531,1532,1533,1534,1535,1536,1537,1538,1539,1540,1541,1542,1543,1544,1545,1546,1547,1548,1549,1550,1551,1552,1553]
        if dictMe[method] == 'b3':
            indexb3 = [1554,1555,1556,1557,1558,1559,1560,1561,1562]
        if dictMe[method] == 'b4':
            indexb4 = [1563,1564]
        if dictMe[method] == 'b5':
            indexb5 = [1565,1566,1567,1568,1569,1570,1571]
        if dictMe[method] == 'b6':
            indexb6 = [1572,1573,1574]
        if dictMe[method] == 'b7':
            indexb7 = list(range(1575,1785,1))
        if dictMe[method] == 'c1':
            indexc1 = list(range(1876,2757,1))
    index = [0] + indexa1 + indexa2 + indexa3 + indexa4 + indexa5 + indexb1 + indexb2 + indexb3 + indexb4 + indexb5 + indexb6 + indexb7 + indexc1
    # dfresult = pd.read_csv(ResfilePath, index_col=0)
    statics.order_DfFile(ResfilePath, InterFilepath).get_interaction_pair()

    dfresult_all = pd.read_csv(ResfilePath)
    dfresult = dfresult_all.iloc[:,index]
    # print('dfresult.head()')
    # print(dfresult.head())
    dfresult.index = dfresult.iloc[:,0].tolist()
    dfresult01 = dfresult.iloc[:,1:]
    # print('dfresult01.head()')
    # print(dfresult01.head())
    dfresult01.columns = list(map(iterlist,index[1:]))
    dfresult02 = dfresult01.dropna(axis='columns', how='any')
    return dfresult02

def iterlist(x):
    y = Dict_methods.Molefeaturenames[x]
    return y
