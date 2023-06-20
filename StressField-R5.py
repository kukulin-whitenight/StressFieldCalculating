import pandas as pd
import openpyxl
import numpy as np
'''
此代码参考文献为：
[1] Li XY, Gao YF, Ge PQ, ZhangL, Bi WB, Meng JF. Nucleation location and propagation direction of radial and median cracks for brittle material in scratching [J]. Ceram. Int., 2019.
[2] Yang X, Gao S. Analysis of the crack propagation mechanism of multiple scratched glass-ceramics by an interference stress field prediction model and experiment [J]. Ceram.Int., 2022
[3] Marshall DB, Lawn BR, Evans AG. Elastic/Plastic indentation damage in ceramics: the lateral crack system [J]. J. Amer. Ceram., 1982
参考文献[1]中，XY平面上的应力通过分别计算Boussinesq、Cerruti和Blister(x<0)的应力，再计算Boussinesq+Cerruti+miu*Blister，其中miu∈[0,1]，0表示材料为脆性去除，1表示材料塑性去除
ZX平面上的应力分别计算Boussinesq、Cerruti和Blister(x>0)的应力,再计算Boussinesq+Cerruti+Blister
计算所需参数：
Sheet='property'  H硬度  v泊松比  E弹性模量
Sheet='condition' a压入深度   fEH选值详见参考文献[1] APPENDIX A     phi磨粒半顶角
                  f the fraction of volume contraction induced by residual stress to the total volume contraction, which is 0 for volume accommodation entirely by densification and 1 for when no densification occors
                  lamda(不用lambda是因为存在代码冲突)无量纲几何系数lamda=2对应于pyramid indenter而lamda=pi对应于conical indenter

注：在塑性变形区h内的应力会远大于其余区域
'''
def read_condition(Filename,Sheetname,Targetname,Mode=0):
    '''
此方法用于读取指定行/列的数据
数据在EXCEL表中的存储方式看我发的模板 ScratchDepthR2-ModuleExcel.xlsx
    :param Filename: 数据文件名
    :param Sheetname: 文件中的工作表名
    :param Targetname: 目标读取的数据索引名
    :param Mode: 读取方式，0：读取行数据；1：读取列数据
    :return:按行/列读取数据，输出对应目标读取的数据索引名的行/列
    '''
    df=pd.read_excel(Filename,sheet_name=Sheetname,index_col=0)
    keys=[]
    values=[]
    if Mode==0:
        keys=df.columns
        values=df.loc[Targetname]
    elif Mode==1:
        keys=df.index
        values=df[Targetname]
    results=dict(zip(keys,values))
    return results

class BoussinesqFieldNormalForce():
    '''
    笛卡儿坐标系，法向力引起的Boussinesq问题解 类
    '''
    def xigma_X(self,P,v,x,y,z,r,rou):
        xigmanx=P / 2 / np.pi * (((1 - 2 * v) / r ** 2 * ((1 - z / rou) * (
                    x ** 2 - y ** 2) / r ** 2 + z * y ** 2 / rou ** 3)) - 3 * z * x ** 2 / rou ** 5)
        return xigmanx
    def xigma_Y(self,P,v,x,y,z,r,rou):
        xigmany = P / 2 / np.pi * (((1 - 2 * v) / r ** 2 * ((1 - z / rou) * (
                y ** 2 - x ** 2) / r ** 2 + z * x ** 2 / rou ** 3)) - 3 * z * y ** 2 / rou ** 5)
        return xigmany
    def xigma_Z(self,P,z,rou):
        xigmanz = -3 * P / 2 / np.pi * z ** 3 / rou ** 5
        return xigmanz
    def tao_XY(self,P,v,x,y,z,r,rou):
        taonxy = P / 2 / np.pi * (((1 - 2 * v) / r ** 2 * ((1 - z / rou) * (2*
                    y * x) / r ** 2 + x * y * z / rou ** 3)) - 3 * x * y * z / rou ** 5)
        return taonxy
    def tao_YZ(self,P,y,z,rou):
        taonyz = -3 * P / 2 / np.pi * y * z ** 2 / rou ** 5
        return taonyz
    def tao_ZX(self,P,x,z,rou):
        taonzx = -3 * P / 2 / np.pi * x * z ** 2 / rou ** 5
        return taonzx
    @classmethod
    def calculating(cls,Property,Condition,Position,Mode):
        '''
        :param Property: 材料性能 H v E  数据类型：字典
        :param Condition: a fEH f phi lamda 划擦条件  数据类型：字典
        :param Position: 位置点 数据类型：字典
        :param Mode: 选填  'XY' / 'YZ' / 'ZX'
        :return: 返回需要的两个正应力和一个对应的切应力
        '''
        v=Property['v']
        H=Property['H']
        phi=Condition['phi']/180*np.pi
        a=Condition['a']
        P=1/2*np.pi*(a*np.tan(phi))**2*H
        x=Position['x']
        y=Position['y']
        z=Position['z']
        r=[]
        rou=[]
        for i in range(len(x)):
            r.append(np.sqrt(x[i]**2+y[i]**2))
            rou.append(np.sqrt(x[i]**2+y[i]**2+z[i]**2))
        xigma1=[]
        xigma2=[]
        tao12=[]
        if Mode=='XY':
            for i in range(len(x)):
                xigma1.append(cls.xigma_X(cls,P,v,x[i],y[i],z[i],r[i],rou[i]))
                xigma2.append(cls.xigma_Y(cls,P, v, x[i], y[i], z[i], r[i], rou[i]))
                tao12.append(cls.tao_XY(cls,P,v,x[i],y[i],z[i],r[i],rou[i]))
        elif Mode=='ZX':
            for i in range(len(z)):
                xigma1.append(cls.xigma_Z(cls,P, z[i],rou[i]))
                xigma2.append(cls.xigma_X(cls,P, v, x[i], y[i], z[i], r[i], rou[i]))
                tao12.append(cls.tao_ZX(cls,P,x[i],z[i],rou[i]))
        return xigma1,xigma2,tao12

    @classmethod
    def calculatingXigma(cls,Property,Condition,Position):
        v=Property['v']
        H=Property['H']
        phi=Condition['phi']/180*np.pi
        a=Condition['a']
        P=1/2*np.pi*(a*np.tan(phi))**2*H
        x=Position['x']
        y=Position['y']
        z=Position['z']
        r=[]
        rou=[]
        for i in range(len(x)):
            r.append(np.sqrt(x[i]**2+y[i]**2))
            rou.append(np.sqrt(x[i]**2+y[i]**2+z[i]**2))
        xigma_x=[]
        xigma_y=[]
        xigma_z=[]
        for i in range(len(x)):
            xigma_x.append(cls.xigma_X(cls,P,v,x[i],y[i],z[i],r[i],rou[i]))
            xigma_y.append(cls.xigma_Y(cls,P, v, x[i], y[i], z[i], r[i], rou[i]))
            xigma_z.append(cls.xigma_Z(cls,P, z[i],rou[i]))
        return xigma_x,xigma_y,xigma_z

class BoussinesqFieldTangentialForce():
    '''
    笛卡儿坐标系 切向力引起的Boussinesq问题解类 也称Cerruti问题解
    '''
    def xigma_X(self,Q,v,x,z,rou):
        xigmatx = -Q / 2 / np.pi * (3 * x ** 3 / rou ** 5 - (1 - 2 * v) * (
                    x / rou ** 3 - 3 * x / rou / (
                        rou + z) ** 2 + x ** 3 / rou ** 3 / (
                                rou + z) ** 2 + 2 * x ** 3 / rou ** 2 / (rou + z) ** 3))
        return xigmatx
    def xigma_Y(self,Q,v,x,y,z,rou):
        xigmaty = -Q / 2 / np.pi * (3 * x * y ** 2 / rou ** 5 - (1 - 2 * v) * (
                    x / rou ** 3 - x / rou / (
                        rou + z) ** 2 + x * y ** 2 / rou ** 3 / (
                                rou + z) ** 2 + 2 * x * y ** 2 / rou ** 2 / (
                                rou + z) ** 3))
        return xigmaty
    def xigma_Z(self,Q,x,z,rou):
        xigmatz = -Q / 2 / np.pi * 3 * x * z ** 2 / rou ** 5
        return xigmatz
    def tao_XY(self,Q,v,x,y,z,rou):
        taotxy = -Q / 2 / np.pi * (3 * x**2 * y / rou ** 5 - (1 - 2 * v) * (
                    y / rou / (rou + z) ** 2 - x ** 2 * y / rou ** 3 / (
                        rou + z) ** 2 - 2 * x ** 2 * y / rou ** 2 / (rou + z) ** 3))
        return taotxy
    def tao_YZ(self,Q,x,y,z,rou):
        taotyz = -Q / 2 / np.pi * 3 * x * y * z / rou ** 5
        return taotyz
    def tao_ZX(self,Q,x,z,rou):
        taotzx = -Q / 2 / np.pi * 3 * x ** 2 * z / rou ** 5
        return taotzx
    @classmethod
    def calculating(cls,Property,Condition,Position,Mode):
        '''

        :param Property: 材料性能  数据类型：字典
        :param Condition: 划擦条件  数据类型：字典
        :param Position: 位置点 数据类型：字典
        :param Mode: 选填  'XY' / 'YZ' / 'ZX'
        :return: 返回需要的两个正应力和一个对应的切应力
        '''
        v=Property['v']
        H=Property['H']
        phi=Condition['phi']/180*np.pi
        a=Condition['a']
        Q=2*(1/np.tan(phi))/np.pi*1/2*np.pi*(a*np.tan(phi))**2*H
        x=Position['x']
        y=Position['y']
        z=Position['z']
        rou = []
        for i in range(len(x)):
            rou.append(np.sqrt(x[i] ** 2 + y[i] ** 2 + z[i] ** 2))
        xigma1=[]
        xigma2=[]
        tao12=[]
        if Mode=='XY':
            for i in range(len(x)):
                xigma1.append(cls.xigma_X(cls,Q,v,x[i],z[i],rou[i]))
                xigma2.append(cls.xigma_Y(cls,Q, v, x[i], y[i], z[i], rou[i]))
                tao12.append(cls.tao_XY(cls,Q,v,x[i],y[i],z[i],rou[i]))
        elif Mode=='ZX':
            for i in range(len(z)):
                xigma1.append(cls.xigma_Z(cls,Q,x[i], z[i],rou[i]))
                xigma2.append(cls.xigma_X(cls,Q, v, x[i], z[i], rou[i]))
                tao12.append(cls.tao_ZX(cls,Q,x[i],z[i],rou[i]))
        return xigma1,xigma2,tao12

    @classmethod
    def calculatingXigma(cls,Property,Condition,Position):
        v=Property['v']
        H=Property['H']
        phi=Condition['phi']/180*np.pi
        a=Condition['a']
        Q=2*(1/np.tan(phi))/np.pi*1/2*np.pi*(a*np.tan(phi))**2*H
        x=Position['x']
        y=Position['y']
        z=Position['z']
        rou = []
        for i in range(len(x)):
            rou.append(np.sqrt(x[i] ** 2 + y[i] ** 2 + z[i] ** 2))
        xigma_x=[]
        xigma_y=[]
        xigma_z=[]
        for i in range(len(x)):
            xigma_x.append(cls.xigma_X(cls, Q, v, x[i], z[i], rou[i]))
            xigma_y.append(cls.xigma_Y(cls, Q, v, x[i], y[i], z[i], rou[i]))
            xigma_z.append(cls.xigma_Z(cls, Q, x[i], z[i], rou[i]))
        return xigma_x,xigma_y,xigma_z

class BlisterField():
    '''
    残余应力解类（？）
    '''
    def xigma_X(self,BS,v,x,y,z,rou):
        xigmax = 2 * BS / (y ** 2 + z ** 2) ** 2 * (
                    -2 * v * (y ** 2 - z ** 2) + x / rou ** 5 * (
                        2 * v * x ** 4 * y ** 2 - 2 * x ** 2 * y ** 4 + 6 * v * x ** 2 * y ** 4 - 2 * y ** 6 + 4 * v * y ** 6 - 2 * v * x ** 4 * z ** 2 - 4 * x ** 2 * y ** 2 * z ** 2 + 2 * v * x ** 2 * y ** 2 * z ** 2 - 3 * y ** 4 * z ** 2 + 6 * v * y ** 4 * z ** 2 - 2 * x ** 2 * z ** 4 - 4 * v * x ** 2 * z ** 4 + z ** 6 - 2 * v * z ** 6))
        return xigmax
    def xigma_Y(self,BS,v,x,y,z,rou):
        xigmay = 2 * BS / (y ** 2 + z ** 2) ** 3 * (
                    -2 * y ** 2 * (y ** 2 - 3 * z ** 2) + x / rou ** 5 * (
                        2 * x ** 4 * y ** 4 - 2 * v * x ** 2 * y ** 6 + 6 * x ** 2 * y ** 6 + 4 * y ** 8 - 2 * v * y ** 8 - 6 * x ** 4 * y ** 2 * z ** 2 - 7 * x ** 2 * y ** 4 * z ** 2 - 6 * v * x ** 2 * y ** 4 * z ** 2 - 2 * y ** 6 * z ** 2 - 8 * v * y ** 6 * z ** 2 - 12 * x ** 2 * y ** 2 * z ** 4 - 6 * v * x ** 2 * y ** 2 * z ** 4 - 15 * y ** 4 * z ** 4 - 12 * v * y ** 4 * z ** 4 + x ** 2 * z ** 6 - 2 * v * x ** 2 * z ** 6 - 8 * y ** 2 * z ** 6 - 8 * v * y ** 2 * z ** 6 + z ** 8 - 2 * v * z ** 8))
        return xigmay
    def xigma_Z(self,BS,x,y,z,rou):
        xigmaz = 2 * BS * z ** 2 / (y ** 2 + z ** 2) ** 3 * (
                    2 * (z ** 2 - 3 * y ** 2) + x / rou ** 5 * (
                        6 * x ** 4 * y ** 2 + 15 * x ** 2 * y ** 4 + 9 * y ** 6 - 2 * x ** 4 * z ** 2 + 10 * x ** 2 * y ** 2 * z ** 2 + 12 * y ** 4 * z ** 2 - 5 * x ** 2 * z ** 4 - 3 * y ** 2 * z ** 4 - 6 * z ** 6))
        return xigmaz
    def tao_XY(self,BS,v,x,y,z,rou):
        taoxy = -2 * BS * y / rou ** 5 * (2 * (1 - v) * x ** 2 + 2 * (
                    1 - v) * y ** 2 - z ** 2 - 2 * v * z ** 2)
        return taoxy
    def tao_ZX(self,BS,x,y,z,rou):
        taoxz = -2 ** BS * z / rou ** 5 * (2 * x ** 2 + 2 * y ** 2 - z ** 2)
        return taoxz
    def tao_YZ(self,BS,x,y,z,rou):
        taoyx = 2 * BS * y * z / (y ** 2 + z ** 2) ** 3 * (
                    -4 * (y ** 2 - z ** 2) + x / rou ** 5 * (
                        4 * x ** 4 * y ** 2 + 10 * x ** 2 * y ** 4 + 6 * y ** 6 - 4 * x ** 4 * z ** 2 + 3 * y ** 4 * z ** 2 - 10 * x ** 2 * z ** 4 - 12 * y ** 2 * z ** 4 - 9 * z ** 6))
        return taoyx
#x>0部分的的Blister解
    def xigma_X_R(self,BI,v,x,y,z,r,rou):
        xigmax = 2 * BI * ((1 - 2 * v) * (x ** 2 + 2 * y ** 2) / r ** 2 / rou ** 3 + (
                    -5 + 4 * v) * x ** 2 / rou ** 5 - (1 - 2 * v) * (
                                            2 * x ** 2 + 3 * y ** 2) * z ** 2 / r ** 2 / rou ** 5 + 15 * x ** 2 * z ** 2 / rou ** 7)
        return xigmax

    def xigma_Y_R(self,BI,v,x,y,z,r,rou):
        xigmay = 2 * BI * ((1 - 2 * v) * (2 * x ** 2 + y ** 2) / r ** 2 / rou ** 3 + (
                -5 + 4 * v) * y ** 2 / rou ** 5 - (1 - 2 * v) * (
                                        3 * x ** 2 + 2 * y ** 2) * z ** 2 / r ** 2 / rou ** 5 + 15 * y ** 2 * z ** 2 / rou ** 7)
        return xigmay
    def xigma_Z_R(self,BI,z,rou):
        xigmaz = -2 * BI * 3 * z ** 2 / rou ** 5 * (3 - 5 * z ** 2 / rou ** 2)
        return xigmaz
    def tao_XY_R(self,BI,v,x,y,z,r,rou):
        taoxy = -2 * BI * ((1 - 2 * v) * x * y / r ** 2 / rou ** 3 - (
                    -5 + 4 * v) * x * y / rou ** 5 - (
                                            1 - 2 * v) * x * y * z **2 / r ** 2 / rou ** 5 - 15 * x * y * z ** 2 / rou ** 7)
        return taoxy
    def tao_YZ_R(self,BI,y,z,rou):
        taoxz = -2 * BI * (6 * y * z / rou ** 5 - 15 * y * z ** 3 / rou ** 7)
        return taoxz
    def tao_ZX_R(self,BI,x,z,rou):
        taoyx = -2 * BI * (6 * x * z / rou ** 5 - 15 * x * z ** 3 / rou ** 7)
        return taoyx

    @classmethod
    def calculating(cls,Property,Condition,Position,Mode):
        '''

        :param Property: 材料性能  数据类型：字典
        :param Condition: 划擦条件  数据类型：字典
        :param Position: 位置点 数据类型：字典
        :param Mode: 选填  'XY' / 'YZ' / 'ZX'
        :return: 返回需要的两个正应力和一个对应的切应力
        '''
        H=Property['H']
        v=Property['v']
        E=Property['E']
        phi=Condition['phi']/180*np.pi
        a=Condition['a']
        fEH=Condition['fEH']
        lamda=Condition['lamda']
        f=Condition['f']
        V=2*a**3/3/np.tan(phi) #压痕体积
        Q = 2 * (1 / np.tan(phi)) / np.pi * 1 / 2 * np.pi * (a * np.tan(phi)) ** 2 * H
        BS=fEH * 3 / 4 / np.pi / lamda / (1 - 2 * v) / (1 + v) * 1 / np.tan(phi) * Q
        BI = 3 * E / 4 / np.pi / (1 - 2 * v) / (1 + v) * f * V
        x=Position['x']
        y=Position['y']
        z=Position['z']
        r=[]
        rou=[]
        for i in range(len(x)):
            r.append(np.sqrt(x[i]**2+y[i]**2))
            rou.append(np.sqrt(x[i]**2+y[i]**2+z[i]**2))
        xigma1=[]
        xigma2=[]
        tao12=[]
        if Mode=='XY':
            for i in range(len(x)):
                xigma1.append(cls.xigma_X(cls, BS, v, x[i], y[i], z[i], rou[i]))
                xigma2.append(cls.xigma_Y(cls, BS, v, x[i], y[i], z[i], rou[i]))
                tao12.append(cls.tao_XY(cls, BS, v, x[i], y[i], z[i], rou[i]))
        elif Mode=='ZX':
            for i in range(len(x)):
                xigma1.append(cls.xigma_Z_R(cls, BI, z[i], rou[i]))
                xigma2.append(cls.xigma_X_R(cls, BI, v, x[i], y[i], z[i], r[i],rou[i]))
                tao12.append(cls.tao_ZX_R(cls, BI, x[i], z[i], rou[i]))
        return xigma1,xigma2,tao12

    @classmethod
    def calculatingXigma(cls,Property,Condition,Position):
        H = Property['H']
        v = Property['v']
        E = Property['E']
        phi = Condition['phi'] / 180 * np.pi
        a = Condition['a']
        f = Condition['f']
        V = 2 * a ** 3 / 3 / np.tan(phi)  # 压痕体积
        BI = 3 * E / 4 / np.pi / (1 - 2 * v) / (1 + v) * f * V
        x = Position['x']
        y = Position['y']
        z = Position['z']
        r = []
        rou = []
        for i in range(len(x)):
            r.append(np.sqrt(x[i] ** 2 + y[i] ** 2))
            rou.append(np.sqrt(x[i] ** 2 + y[i] ** 2 + z[i] ** 2))
        xigma_x=[]
        xigma_y=[]
        xigma_z=[]
        for i in range(len(x)):
            xigma_x.append(cls.xigma_X_R(cls,BI,v,x[i],y[i],z[i],r[i],rou[i]))
            xigma_y.append(cls.xigma_Y_R(cls,BI,v,x[i],y[i],z[i],r[i],rou[i]))
            xigma_z.append(cls.xigma_Z_R(cls, BI, z[i], rou[i]))
        return xigma_x,xigma_y,xigma_z

def stressFiledCalculating(Property,Condition,Position,miu,Mode):
    '''
    :param Property: 材料性能  数据类型：字典
    :param Condition: 划擦条件  数据类型：字典
    :param Position: 位置点 数据类型：字典
    :param miu:应力场叠加计算中的Blister场计算系数，miu∈[0,1]，0表示材料为脆性去除，1表示材料塑性去除。
    :param Mode: 选填  'XY' / 'ZX'
    :return: 返回需要的两个正应力和一个对应的切应力  xigma1 xigma2 tao12  计算'XY'平面时xigma1为xigma_X,xigma2为xigma_Y 而计算'ZX'时xigma1为xigma_Z,xigma2为xigma_X
    '''

    BoussinesqN_xigma1, BoussinesqN_xigma2, BoussinesqN_tao12 = BoussinesqFieldNormalForce.calculating(Property,
                                                                                                       Condition,
                                                                                                       Position, Mode)
    BoussinesqT_xigma1, BoussinesqT_xigma2, BoussinesqT_tao12 = BoussinesqFieldTangentialForce.calculating(Property,
                                                                                                           Condition,
                                                                                                           Position,
                                                                                                           Mode)
    Blister_xigma1, Blister_xigma2, Blister_tao12 = BlisterField.calculating(Property, Condition, Position, Mode)
    xigma1 = []
    xigma2 = []
    tao12 = []
    if Mode=='XY':
        for i in range(len(BoussinesqN_xigma1)):
            xigma1.append(BoussinesqN_xigma1[i] + BoussinesqT_xigma1[i] + miu*Blister_xigma1[i])
            xigma2.append(BoussinesqN_xigma2[i] + BoussinesqT_xigma2[i] + miu*Blister_xigma2[i])
            tao12.append(BoussinesqN_tao12[i] + BoussinesqT_tao12[i] + miu*Blister_tao12[i])
    elif Mode=='ZX':
        miu=1
        for i in range(len(BoussinesqN_xigma1)):
            xigma1.append(BoussinesqN_xigma1[i] + BoussinesqT_xigma1[i] + miu*Blister_xigma1[i])
            xigma2.append(BoussinesqN_xigma2[i] + BoussinesqT_xigma2[i] + miu*Blister_xigma2[i])
            tao12.append(BoussinesqN_tao12[i] + BoussinesqT_tao12[i] + miu*Blister_tao12[i])
    return xigma1,xigma2,tao12

def dataProcessing(Property,Condition,Position,miu,Mode='XY'):
    '''
    坐标点通过塑性变形区域半径归一化，计算最大主应力并归一化
    :param Property: 材料性能  数据类型：字典
    :param Condition: 划擦条件  数据类型：字典
    :param Position: 位置点 数据类型：字典
    :param miu:应力场叠加计算中的Blister场计算系数，miu∈[0,1]，0表示材料为脆性去除，1表示材料塑性去除。
    :param Mode: 选填  'XY' / 'ZX'
    :return: 返回需要的两个正应力和一个对应的切应力
    '''
    x=Position['x']
    y=Position['y']
    z=Position['z']
    phi=Condition['phi']/180*np.pi
    E=Property['E']
    H=Property['H']
    a=Condition['a']
    P = 1 / 2 * np.pi * (a * np.tan(phi)) ** 2 * H
    h=0.226 * (1 / np.tan(phi)) ** (1 / 3) * E ** (3 / 4) / H * P ** (1 / 2)  #塑性变形区半径 参考文献[3]
    print('P:{}'.format(P))
    print('h:{}'.format(h))
    x_Normalized=[var/h for var in x]
    y_Normalized=[var/h for var in y]
    z_Normalized=[var/h for var in z]
    xigma_1,xigma_2,tao_12=stressFiledCalculating(Property,Condition,Position,miu,Mode)
    xigma_MAX=[]
    print('calculation star')
    for i in range(len(xigma_1)):
        xigma_MAX.append(((xigma_1[i]+xigma_2[i])/2+np.sqrt(((xigma_1[i]-xigma_2[i])/2)**2+tao_12[i]**2)))
    xigma_MAX_Normalized=[2*np.pi*h**2*var/P for var in xigma_MAX]
    ResultDict={}
    ResultDict['x']=x
    ResultDict['y'] = y
    ResultDict['z'] = z
    ResultDict['xigma_1'] = xigma_1
    ResultDict['xigma_2']=xigma_2
    ResultDict['tao_12']=tao_12
    ResultDict['xigma_MAX']=xigma_MAX
    ResultDict['x_Normalized']=x_Normalized
    ResultDict['y_Normalized']=y_Normalized
    ResultDict['z_Normalized']=z_Normalized
    ResultDict['xigma_MAX_Normalized']=xigma_MAX_Normalized
    print('calculation end')
    return ResultDict

def stressFieldcalculatingXigma(Property,Condition,Position):
    BoussinesqN_xigma_x, BoussinesqN_xigma_y, BoussinesqN_xigma_z = BoussinesqFieldNormalForce.calculatingXigma(Property,Condition,Position)
    BoussinesqT_xigma_x, BoussinesqT_xigma_y, BoussinesqT_xigma_z = BoussinesqFieldTangentialForce.calculatingXigma(Property,Condition,Position)
    Blister_xigma_x, Blister_xigma_y, Blister_xigma_z = BlisterField.calculatingXigma(Property, Condition, Position)
    xigma_x=[]
    xigma_y=[]
    xigma_z=[]
    for i in range(len(BoussinesqN_xigma_x)):
        xigma_x.append(BoussinesqN_xigma_x[i]+BoussinesqT_xigma_x[i]+Blister_xigma_x[i])
        xigma_y.append(BoussinesqN_xigma_y[i] + BoussinesqT_xigma_y[i] + Blister_xigma_y[i])
        xigma_z.append(BoussinesqN_xigma_z[i] + BoussinesqT_xigma_z[i] + Blister_xigma_z[i])
    return xigma_x,xigma_y,xigma_z

def dataProcessingXigma(Property,Condition,Position):
    '''
    坐标点通过塑性变形区域半径归一化，计算最大主应力并归一化
    :param Property: 材料性能  数据类型：字典
    :param Condition: 划擦条件  数据类型：字典
    :param Position: 位置点 数据类型：字典
    :param miu:应力场叠加计算中的Blister场计算系数，miu∈[0,1]，0表示材料为脆性去除，1表示材料塑性去除。
    :param Mode: 选填  'XY' / 'ZX'
    :return: 返回需要的两个正应力和一个对应的切应力
    '''
    x=Position['x']
    y=Position['y']
    z=Position['z']
    phi=Condition['phi']/180*np.pi
    E=Property['E']
    H=Property['H']
    a=Condition['a']
    P = 1 / 2 * np.pi * (a * np.tan(phi)) ** 2 * H
    h=0.226 * (1 / np.tan(phi)) ** (1 / 3) * E ** (3 / 4) / H * P ** (1 / 2)  #塑性变形区半径 参考文献[3]
    print('P:{}'.format(P))
    print('h:{}'.format(h))
    x_Normalized=[var/h for var in x]
    y_Normalized=[var/h for var in y]
    z_Normalized=[var/h for var in z]
    xigma_x,xigma_y,xigma_z=stressFieldcalculatingXigma(Property,Condition,Position)
    print('calculation star')
    xigma_x_Normalized=[2*np.pi*h**2*var/P for var in xigma_x]
    xigma_y_Normalized = [2 * np.pi * h ** 2 * var / P for var in xigma_y]
    xigma_z_Normalized = [2 * np.pi * h ** 2 * var / P for var in xigma_z]
    ResultDict={}
    ResultDict['x']=x
    ResultDict['y'] = y
    ResultDict['z'] = z
    ResultDict['xigma_x'] = xigma_x
    ResultDict['xigma_y']=xigma_y
    ResultDict['xigma_z']=xigma_z
    ResultDict['x_Normalized']=x_Normalized
    ResultDict['y_Normalized']=y_Normalized
    ResultDict['z_Normalized']=z_Normalized
    ResultDict['xigma_x_Normalized'] = xigma_x_Normalized
    ResultDict['xigma_y_Normalized'] = xigma_y_Normalized
    ResultDict['xigma_z_Normalized'] = xigma_z_Normalized
    print('calculation end')
    return ResultDict

def positionGenerator(arange1,arange2,h,Mode):
    '''
    计算某一参考面上的坐标点，如Mode=='XY',即生成Z=h的XY平面坐标点，默认生成保留小数点后1位
    为了避免生成的坐标点存在应力计算中的奇点，第一个轴生成的点中有0而第二个轴生成的点中无0
    :param arange1: 第一个轴的坐标点范围与生成间隔  [min,max,interval]
    :param arange2: 第二个轴的坐标点范围与生成间隔  [min,max,interval]
    :param h: 第三个轴的坐标点，恒定值
    :param Mode:选填  'XY' / 'YZ' / 'ZX'
    :return: 各坐标轴点 数据类型：字典
    '''
    Position=[]
    a=[round(var,1) for var in np.arange(arange1[0],arange1[1],arange1[2])]  #修改round(b,1)中的1可以改变保留的小数点数
    b=[round(var,1) for var in np.arange(arange2[0],arange2[1],arange2[2])]  #修改round(b,1)中的1可以改变保留的小数点数
    a.remove(0.0)
    b.remove(0.0)
    for i in a:
        for j in b:
            temporary=[]
            temporary.append(round(i,1))#修改round(i,1)中的1可以改变保留的小数点数
            temporary.append(j)
            temporary.append(h)
            Position.append(temporary)
    if [0,0,0] in Position:
        Position.remove([0,0,0])
    tem1=[]
    tem2=[]
    tem3=[]
    for i in range(len(Position)):
        tem1.append(Position[i][0])
        tem2.append(Position[i][1])
        tem3.append(Position[i][2])
    position_point=[]
    position_point.append(tem1)
    position_point.append(tem2)
    position_point.append(tem3)
    keys=[]
    if Mode=='XY':
        keys = ['x', 'y', 'z']
    elif Mode=='YZ':
        keys=['z','y','x']
    elif Mode=='ZX':
        keys=['x','z','y']
    position_dict=dict(zip(keys,position_point))
    print('Position point generated')
    return position_dict

def report_data(Filename,Sheetname,Datadict):
    '''
    将数据输出至EXCEL
    :param Filename:
    :param Sheetname:
    :param Datadict:
    :return:
    '''
    df=pd.DataFrame()
    for k in Datadict.keys():
        df=pd.concat([df,pd.DataFrame({k:Datadict[k]})],axis=1)
    with pd.ExcelWriter(Filename,engine='openpyxl',mode='a') as write:
        df.to_excel(write,sheet_name=Sheetname,index=False)
if __name__=='main':
    #数据存储位置
    Filename=r'C:\Users\Desktop\StressField.xlsx'
    #读取数据  性能、划擦条件
    material_property=read_condition(Filename,'property','GaAs')
    scratch_condition=read_condition(Filename,'condition','set1')
    #生产若干个Z=0，XY平面上的坐标点
    position_point1=positionGenerator([-5,5,0.1],[-5,5,0.1],0,'XY')
    #计算XY平面应力场结果
    ResultsDict1=dataProcessing(material_property,scratch_condition,position_point1,0,Mode='XY')
    #输出结果  文献[1]中 绘图选择  x_normalized(x/h)、y_normalized(y/h)、xigma_max_normalized(2*np.pi*h**2*xigma_max/P)
    report_data(Filename,'Results_set1_GaAs_XY',ResultsDict1)

    #生产若干个Y=0，ZX平面上的坐标点
    position_point2=positionGenerator([-5,5,0.1],[-5,5,0.1],0,'ZX')
    #计算ZX应力场结果
    ResultsDict2=dataProcessingXigma(material_property,scratch_condition,position_point2)
    #输出结果  文献[1]中 绘图选择  x_normalized(x/h)、z_normalized(z/h)、xigma_y_normalized(2*np.pi*h**2*xigma_y/P)
    report_data(Filename,'Results_set1_GaAs_ZX',ResultsDict2)
