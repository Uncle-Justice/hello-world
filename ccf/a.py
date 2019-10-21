# 导入工具包
import numpy as np
import matplotlib.pyplot as plt
import argparse
import cv2
import scipy.signal
# from PIL import Image

# 显示时间
debugWP = True
debugID = False
show = False
showID = False
showprint = False
drawHist = False
getWP =False
time = 0
wpBgRatio = 0.21 #水印最低灰度/背景灰度 比值
wpGrad = 0.147 #水印边缘灰度梯度 
wpTxDt = 0.08 #水印与被覆盖字体灰度偏差/背景
exposure = 0.8  #曝光：range from 1 to 2
bgx , bgy = 6 , 6   #背景色取点
stdArea = 126000    #身份证区域标准面积
waterArea = 8300    #水印标准面积
areaOffSet = 0.08   #面积偏差容忍度
wareaOffSet = 0.2
wSimilar = 0.22
sharpen = 1   #锐化程度: range from 1 to 2
wpCntIn = np.array([[171, 46], [  0, 46], [  0,  0], [171 , 0]],np.int64)
wpCntOut = np.array([[  0,  0], [179 , 0],[179, 51], [  0, 51]],np.float32)
path = 'C:/Users/SoWhat/Desktop/CCF/'
err = 0

def findMedian(self,lst):
        #先将列表进行排序
        lst.sort()
        half = len(lst)//2
        #得到中间序列，~half为负索引，列表元素可能为偶数，需要获取中间两个数
        #转化成float，中位数可能为浮点数类型，如测试用例
        median = (float(lst[half])+float(lst[~half]))/2
        return median

def order_points(pts):
    # 一共4个坐标点
    rect = np.zeros((4, 2), dtype="float32")

    # 按顺序找到对应坐标0123分别是 左上，右上，右下，左下
    # 计算左上，右下
    s = pts.sum(axis=1)
    rect[0] = pts[np.argmin(s)]
    rect[2] = pts[np.argmax(s)]

    # 计算右上和左下
    diff = np.diff(pts, axis=1)
    rect[1] = pts[np.argmin(diff)]
    rect[3] = pts[np.argmax(diff)]

    return rect

def four_point_transform(image, pts, mode = 0 , w = 0, h = 0):
    #mode == 0 从image里截出pts并仿射至矩形
    #mode == 1 把image仿射至新图中的pts位置并设置长宽
    # 获取输入坐标点
    rect = order_points(pts)
    (tl, tr, br, bl) = rect
    if mode == 0:

        # 计算输入的w和h值
        widthA = np.sqrt(((br[0] - bl[0]) ** 2) + ((br[1] - bl[1]) ** 2))
        widthB = np.sqrt(((tr[0] - tl[0]) ** 2) + ((tr[1] - tl[1]) ** 2))
        maxWidth = max(int(widthA), int(widthB))

        heightA = np.sqrt(((tr[0] - br[0]) ** 2) + ((tr[1] - br[1]) ** 2))
        heightB = np.sqrt(((tl[0] - bl[0]) ** 2) + ((tl[1] - bl[1]) ** 2))
        maxHeight = max(int(heightA), int(heightB))

        # 变换后对应坐标位置
        dst = np.array([
            [0, 0],
            [maxWidth - 1, 0],
            [maxWidth - 1, maxHeight - 1],
            [0, maxHeight - 1]], dtype="float32")

        # 计算变换矩阵
        M = cv2.getPerspectiveTransform(rect, dst)
        warped = cv2.warpPerspective(image, M, (maxWidth, maxHeight))
        # 返回变换后结果
        return warped
    elif mode == 1:
        # 计算变换矩阵
        M = cv2.getPerspectiveTransform(wpCntOut, order_points(pts))
        warped = cv2.warpPerspective(image, M, (w, h))
        return warped

    else:
        return

def resize(image, width=None, height=None, inter=cv2.INTER_AREA):
    dim = None
    (h, w) = image.shape[:2]
    if width is None and height is None:
        return image
    if width is None:
        r = height / float(h)
        dim = (int(w * r), height)
    elif height is None:
        r = width / float(w)
        dim = (width, int(h * r))
    else:
        dim = (width, height)
    resized = cv2.resize(image, dim, interpolation=inter)
    return resized

def calHistogram(img):
    if(len(img.shape) != 2):
        print("img size error")
        return None
    histogram = {}
    for i1 in range(img.shape[0]):
        for i2 in range(img.shape[1]):
            if histogram.get(img[i1][i2]) is None:
                histogram[img[i1][i2]] = 0
            histogram[img[i1][i2]] += 1
    # normalize
    for key in histogram:
        histogram[key] = float(histogram[key]) * 100 / (img.shape[0]*img.shape[1])
    return histogram

def drawHistoGram(img):
    histogram = calHistogram(img)
    plt.figure()
    #设置x轴的最小值，最大值。y轴的最小值，最大值
    plt.axis([0, 256, 0, 10])
    #显示网格线
    plt.grid(True)
    #key正好就是灰度
    keys = histogram.keys()
    #value是灰度的像素数量，这里是归一化之后的
    values = histogram.values()
    #这里正式绘制直方图
    plt.bar(tuple(keys), tuple(values))
    plt.show()

def findWP(a1,m1,pts,waterCnt):
    #按面积确定水印轮廓
    if  abs( a1 - waterArea ) < wareaOffSet * waterArea \
        and cv2.matchShapes(wpCntIn, pts, 1, 0.0) < wSimilar:
        #去重
        if len(waterCnt) >= 1:
            lpts = waterCnt.pop()
            m2 = cv2.moments(lpts)
            a2 = cv2.contourArea(lpts)

            #重心远离：判断为另一轮廓
            if abs(m1['m01']/m1['m00']-m2['m01']/m2['m00']) > 100 \
                and len(waterCnt) == 0:
                waterCnt.append(lpts)
                waterCnt.append(pts)
                if debugWP:
                    print("select" + str(a1) +" and " + str(a2))
            else:
                waterCnt.append(lpts)
        else:
            waterCnt.append(pts)
            if debugWP:
                print("First: " + str(a1))

def findContWithArea(img,cnts,stdArea,areaOffSet,CntList,debug = False):
    for pts in cnts:
        # 计算面积及重心

        m1 = cv2.moments(pts)
        a1 = cv2.contourArea(pts)

        if debug:
            print(a1)
            copy = img.copy()
            cv2.drawContours(copy, [pts], 0, (0, 255, 0), 2)
            cv2.imshow("IDcardOutline", copy)  
            cv2.waitKey(time)
            cv2.destroyAllWindows()

        #按面积确定轮廓
        if  abs( a1 - stdArea ) < areaOffSet * stdArea :
            #去重
            if len(CntList) >= 1:
                #print('m1 ' + str(m1['m01']/m1['m00']))
                lpts = CntList.pop()
                m2 = cv2.moments(lpts)
                a2 = cv2.contourArea(lpts)
                #print('m2 ' + str(m2['m01']/m2['m00']))

                #重心远离：判断为另一面轮廓
                if abs(m1['m01']/m1['m00']-m2['m01']/m2['m00']) > 150 \
                    and len(CntList) == 0:
                    CntList.append(lpts)
                    CntList.append(pts)
                else:
                    CntList.append(lpts)
            else:
                CntList.append(pts)

def rectBound(cnts,CntList,rev = False):
    for c in cnts:
        rect = cv2.minAreaRect(c)
        box = cv2.boxPoints(rect)
        pts = np.int0(box)
        CntList.append(pts)
    return sorted(CntList, key=cv2.contourArea, reverse = rev)

def getBgG(img):
    alpha = 0
    Hist = cv2.calcHist([img], [0], None, [256], [0,256])
    for i in range(1,256):
        m = max(Hist)
        if Hist[i] == m:
            alpha = i
            print("alpha: " + str(alpha))
            break
    return alpha

def cutHist(img, min , max , delta, mode = 0, mask = None):
    #TODO:thresh - thresh bitwise 25 -106
    #将min -max 拉高delta
    ret1,threshMin = cv2.threshold(img, min, 255, cv2.THRESH_BINARY)
    ret2,threshMax = cv2.threshold(img, max, 255, cv2.THRESH_BINARY_INV)
    part = cv2.bitwise_and(threshMax,threshMin,mask = mask)
    if mode == 0:
        res = cv2.add(img, part // 255 * delta)
    #将min -max 拉到delta
    if mode == 1:
        cv2.bitwise_xor(img.copy(), img.copy(),img, mask = part) #img的part 部分变为0
        #cv2.imshow("res", img)
        #cv2.waitKey(time)
        #cv2.destroyAllWindows()
        res = cv2.add(img, part // 255 * delta)
    return res

def transform(imagePath,show=False):
    # 读取输入
    image = cv2.imread(imagePath)

    ###### 预处理 ##########
    #灰度
    gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
    #clahe = cv2.createCLAHE(clipLimit=3.0, tileGridSize=(5,5)) 
    #res_clahe = clahe.apply(gray.copy())
    # gray = cv2.GaussianBlur(gray, (5, 5), 0)
    #gray = cutHist(gray,gray[6][6] - 100, gray[6][6] - 55, 80)

    if drawHist:
        cv2.imshow("gray", gray)
        cv2.waitKey(time)
        cv2.destroyAllWindows()
        drawHistoGram(gray)
        return

    #锐化：也叫边缘增强
    kernel = np.array([[-0.7 * sharpen, -sharpen, -0.7 * sharpen], [-sharpen, 7 * sharpen  , -sharpen], [-0.7 * sharpen, -sharpen, -0.7 * sharpen]], np.float32)
    #np.array([[0, -1.4, 0], [-1.4, 6.2, -1.4], [0, -1.4, 0]], np.float32) 
    dst = cv2.filter2D(gray, -1, kernel=kernel)

    
    bg = dst[bgx][bgy]
    
    #二值化
    ret, thresh = cv2.threshold(dst, bg, 255, cv2.THRESH_BINARY)
    
    if debugID:
        cv2.imshow("threshID", thresh)
        cv2.waitKey(time)
        cv2.destroyAllWindows()

    #水印预处理 64 - 145(bg170, peak 164)
    wp = cv2.imread(path + 'waterprint/thresh2.jpg',cv2.IMREAD_GRAYSCALE)
    r, wp = cv2.threshold(wp, 127, 255, cv2.THRESH_BINARY_INV)
    wp = wp // 255
    if getWP:
        retW, wp = cv2.threshold(wp, 178, 255, cv2.THRESH_BINARY)
        cv2.imshow("wp",wp)
        cv2.imwrite(path + 'waterprint/thresh2.jpg',wp)
        wContour = cv2.findContours(wp, cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)[1]
        wContour = sorted(wContour, key=cv2.contourArea, reverse=True)[:2]
        rect = cv2.minAreaRect(wContour[0])
        box = cv2.boxPoints(rect)
        pts = np.int0(box)
        print(pts)
        cv2.waitKey(time)

    
    #clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(4,4)) 
    #res_clahe = clahe.apply(gray.copy())
    
    
    
    ###边缘（未使用）
    #edged = cv2.Canny(thresh, 120, 170)
    
    # 展示预处理结果
    if show:
        print("STEP 1: 边缘检测")
        result = np.hstack((dst,thresh))
        cv2.imshow("Result", result)
        # cv2.imshow("Edged", edged)
        cv2.waitKey(time)
        cv2.destroyAllWindows()

    ## card轮廓检测 #####################################################################################################################
    #这里获取所有轮廓的坐标集合
    cnts = cv2.findContours(thresh.copy(), cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)[1]
    cnts = sorted(cnts, key=cv2.contourArea, reverse=True)[:12]
    
    screenCnt = []
    # 遍历轮廓，并排序
    # 计算IDcard轮廓矩形近似
    cnts = rectBound(cnts,screenCnt,True)
    screenCnt.clear()
    #按面积和重心确定两个轮廓
    findContWithArea(image,cnts,stdArea,areaOffSet,screenCnt,debugID)
    

    ### 水印轮廓 #########################################################################################################################
    
    wpCnt = []
    for sc in screenCnt:
        #sc = np.array(sc)
        waterCnt = []
        wraped = four_point_transform(gray, sc)
        #print(wraped.shape)
        #wraped = np.rollaxis(wraped,1)
        #print(wraped.shape)

        #十字方向边缘加强
        clahe = cv2.createCLAHE(clipLimit = 3.0, tileGridSize=(5,5)) 
        res_clahe = clahe.apply(wraped.copy())
        threshR = cv2.filter2D(res_clahe, -1, kernel = np.array([[0, -1.3, 0], [-1.3, 6, -1.3], [0, -1.3, 0]], np.float32))
        if debugWP:    
            cv2.imshow("sharpen", threshR)
            cv2.waitKey(time)
            cv2.destroyAllWindows()
        #二值化
        #retR, threshR = cv2.threshold(threshR, 20 , 255, cv2.THRESH_TOZERO)
        kernel_hat = np.ones((3,3),np.uint8)
        alpha = getBgG(threshR)
            
        retR, threshR = cv2.threshold(threshR, alpha - 50, 255, cv2.THRESH_BINARY)
        #threshR = cv2.dilate(threshR,kernel_hat,iterations = 1)
        #threshR = cv2.erode(threshR,kernel_hat,iterations = 1)
        
        if debugWP:    
            cv2.imshow("binary", threshR)
            cv2.waitKey(time)
            cv2.destroyAllWindows()

        #形态学膨胀与梯度操作
        kernel_v = np.array([[0,1,0],[0,1,0],[0,1,0]],np.uint8)  
        kernel_h = np.array([[0,0,0],[1,1,1],[0,0,0]],np.uint8) 
        threshR = cv2.dilate(threshR,kernel_hat,iterations = 1)
        
        cntsR = cv2.findContours(threshR.copy(), cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)[1]
        cntsR = sorted(cntsR, key=cv2.contourArea, reverse=True)[:8]
    
        # 计算waterprint轮廓矩形近似
        cntsR = rectBound(cntsR,waterCnt)
        waterCnt.clear()
        
        for pts in cntsR:
            # 计算面积及重心
            m1 = cv2.moments(pts)
            a1 = cv2.contourArea(pts)

            if debugWP:
                print(a1)
                print(cv2.matchShapes(wpCntIn, pts, 1, 0.0))
                copy = wraped.copy()
                #cv2.drawContours(copy, [wpCntIn], 0, (0, 200, 255), 2)
                #cv2.drawContours(copy, [pts], 0, (0, 255, 0), 2)
                wpDebug = four_point_transform(copy, pts)
                cv2.imshow("WaterPrintOutline", wpDebug)  
                cv2.waitKey(time)
                cv2.destroyAllWindows()

            #按面积确定水印轮廓
            findWP(a1,m1,pts,waterCnt)
        
        #备选
        if len(waterCnt) < 1 :
            threshR2 = cv2.morphologyEx(threshR, cv2.MORPH_GRADIENT, kernel_hat)
            if debugWP:    
                cv2.imshow("erode_dilate", threshR2)
                cv2.waitKey(time)
                cv2.destroyAllWindows()
            #threshR = cv2.dilate(threshR,kernel_h,iterations = 3)
            #threshR = cv2.erode(threshR,kernel_h,iterations = 2)

            cntsR = cv2.findContours(threshR2.copy(), cv2.RETR_TREE, cv2.CHAIN_APPROX_SIMPLE)[1]
            cntsR = sorted(cntsR, key=cv2.contourArea, reverse=True)[:8]
            cntsR = rectBound(cntsR,waterCnt)
            waterCnt.clear()
            for pts in cntsR:
                m1 = cv2.moments(pts)
                a1 = cv2.contourArea(pts)
                if debugWP:
                    print(a1)
                    print(cv2.matchShapes(wpCntIn, pts, 1, 0.0))
                    copy = wraped.copy()
                    #cv2.drawContours(copy, [wpCntIn], 0, (0, 200, 255), 2)
                    #cv2.drawContours(copy, [pts], 0, (0, 255, 0), 2)
                    wpDebug = four_point_transform(copy, pts)
                    cv2.imshow("WaterPrintOutline2", wpDebug)  
                    cv2.waitKey(time)
                    cv2.destroyAllWindows()
                findWP(a1,m1,pts,waterCnt)
        #TODO:确定水印外沿
        if len(waterCnt) > 0:
            alpha = getBgG(wraped)
            #输出预处理
            #获取背景灰度与255差值（比）
            bgc = (255 - alpha) / 255
            #np.ones(gray.shape,np.uint8)
            #曝光
            wraped = cv2.add(np.array(wraped * bgc * exposure, np.uint8), wraped)
            #afterProcess = cv2.addWeighted(gray, 1, gray * bgc , exposure / 255, 0)
            #cv2.imshow("曝光", waraped)
            #cv2.waitKey(time)
            #锐化
            wraped = cv2.filter2D(wraped, -1, kernel=np.array([[-0.7 * sharpen, -sharpen, -0.7 * sharpen], [-sharpen, 7.9 * sharpen  , -sharpen], [-0.7 * sharpen, -sharpen, -0.7 * sharpen]], np.float32))
            #cv2.imshow("锐化", wraped)
            #cv2.waitKey(time)
            #r, threshWP = cv2.threshold(wraped,alpha - 1, 255, cv2.THRESH_BINARY)
            if debugWP:    
                cv2.imshow("wraped", wraped)
                cv2.waitKey(time)
                cv2.destroyAllWindows()
            
            alpha = getBgG(wraped)
            ###从里到外四个点，再从外到里四个点，mask一下
            rect = order_points(waterCnt[0])
            (tl, tr, br, bl) = rect
            (TLO, TRO, BRO, BLO) = rect
            (TLI, TRI, BRI, BLI) = rect
            for i in range(-5, 5):
                TLO = tl + np.array([ i,  i])
                #外灰度减内灰度
                deltaTL = int(wraped[int(TLO[1])][int(TLO[0])]) - int(wraped[int(TLO[1]) + 1][int(TLO[0]) + 1])
                if deltaTL > wpGrad * alpha or deltaTL < - wpTxDt * alpha and int(wraped[int(TLO[1])][int(TLO[0])]) > wpBgRatio * alpha:
                    TLO = TLO + np.array([ -1, -1])
                    break
                elif i == 4:
                    TLO = tl + np.array([- 3, - 3])
            for i in range(-5, 5):
                TRO = tr + np.array([  -i,  i])
                deltaTR = int(wraped[int(TRO[1])][int(TRO[0])]) - int(wraped[int(TRO[1]) - 1][int(TRO[0]) + 1])
                if deltaTR > wpGrad * alpha or deltaTR < - wpTxDt * alpha and int(wraped[int(TRO[1])][int(TRO[0])]) > wpBgRatio * alpha:
                    TRO = TRO + np.array([ 1, -1])
                    break
                elif i == 4:
                    TRO = tr + np.array([ 3, - 3])
            for i in range(-5, 5):
                BRO = br + np.array([ -i,  -i])
                deltaBR = int(wraped[int(BRO[1])][int(BRO[0])]) - int(wraped[int(BRO[1]) - 1][int(BRO[0]) - 1])
                if deltaBR > wpGrad * alpha or deltaBR < - wpTxDt * alpha and int(wraped[int(BRO[1])][int(BRO[0]) ]) > wpBgRatio * alpha:
                    BRO = BRO + np.array([ 1,  1])
                    break
                elif i == 4:
                    BRO = br + np.array([ 3,  3])
            for i in range(-5, 5):
                BLO = bl + np.array([ i, - i])
                deltaBL = int(wraped[int(BLO[1])][int(BLO[0])]) - int(wraped[int(BLO[1]) + 1][int(BLO[0]) - 1])
                if deltaBL > wpGrad * alpha or deltaBL < - wpTxDt * alpha and int(wraped[int(BLO[1])][int(BLO[0]) ]) > wpBgRatio * alpha:
                    BLO = BLO + np.array([ -1,  1])
                    break
                elif i == 4:
                    BLO = bl + np.array([ -3, 3])


            for i in range(-5, 5):
                TLI = tl + np.array([ -i,  -i])
                #内灰度减外灰度
                deltaTL = int(wraped[int(TLI[1])][int(TLI[0])]) - int(wraped[int(TLI[1]) - 1][int(TLI[0]) - 1])
                if deltaTL > wpGrad * alpha or deltaTL < -wpTxDt * alpha and int(wraped[int(TLI[1])][int(TLI[0])]) > wpBgRatio * alpha:
                    #TLI = TLI + np.array([ 1, 1])
                    break
                elif i == 4:
                    TLI = tl + np.array([1, 1])
            for i in range(-5, 5):
                TRI = tr + np.array([ i,  -i])
                deltaTR = int(wraped[int(TRI[1])][int(TRI[0])]) - int(wraped[int(TRI[1]) + 1][int(TRI[0]) - 1])
                if deltaTR > wpGrad * alpha or deltaTR < -wpTxDt * alpha and int(wraped[int(TRI[1])][int(TRI[0])]) > wpBgRatio * alpha:
                    #TRI = TRI + np.array([ -1, 1])
                    break
                elif i == 4:
                    TRI = tr + np.array([-1, 1])
            for i in range(-5, 5):
                BRI = br + np.array([ i,  i])
                #内灰度减外灰度
                deltaBR = int(wraped[int(BRI[1])][int(BRI[0])]) - int(wraped[int(BRI[1]) + 1][int(BRI[0]) + 1])
                if deltaBR > wpGrad * alpha or deltaBR < -wpTxDt * alpha and int(wraped[int(BRI[1])][int(BRI[0])]) > wpBgRatio * alpha:
                    #BRI = BRI + np.array([ -1, -1])
                    break
                elif i == 4:
                    BRI = br + np.array([-1, -1])
            for i in range(-5, 5):
                BLI = bl + np.array([ -i,  i])
                #内灰度减外灰度
                deltaBL = int(wraped[int(BLI[1])][int(BLI[0])]) - int(wraped[int(BLI[1]) - 1][int(BLI[0]) + 1])
                if deltaBL > wpGrad * alpha or deltaBL < -wpTxDt * alpha and int(wraped[int(BLI[1])][int(BLI[0])]) > wpBgRatio * alpha:
                    #BLI = BLI + np.array([ 1, -1])
                    break
                elif i == 4:
                    BLI = bl + np.array([1, -1])
            #if abs(int(wraped[int(BLO[1])][int(BLO[0])]) - int(wraped[int(bl[1])][int(bl[0])])) < 0.4 * wraped[int(bl[1])][int(bl[0])]:
            cntO = np.array([TLO, TRO, BRO, BLO],np.int32)
            cntI = np.array([TLI, TRI, BRI, BLI],np.int32)
            mask = np.zeros(wraped.shape,np.uint8)
            maskIn = cv2.bitwise_not(cv2.fillConvexPoly(mask.copy(), cntI, 255))
            maskOut = cv2.fillConvexPoly(mask.copy(), cntO,  255)
            cv2.bitwise_and(maskIn,maskOut.copy(),mask)
            cv2.imshow("mask", cv2.bitwise_and(wraped,wraped,mask = mask))  
            cv2.waitKey(time)
            cv2.destroyAllWindows()
            #cv2.imshow("maskOut", maskOut)  
            #cv2.waitKey(time)
            #cv2.destroyAllWindows()
            if debugWP:
                print("tl " + str(wraped[int(TLO[1])][int(TLO[0])]) \
                + " tr " + str(wraped[int(TRO[1])][int(TRO[0])]) \
                + " br " + str(wraped[int(BRO[1])][int(BRO[0])]) \
                + " bl " + str(wraped[int(BLO[1])][int(BLO[0])]))
                copy = wraped.copy()
                #wpDebug = cv2.circle(copy,(int(TLO[0]),int(TLO[1])),1, color = (0,0,0))
                #wpDebug = cv2.circle(wpDebug,(int(TRO[0]),int(TRO[1])),1, color = (0,0,0))
                #wpDebug = cv2.circle(wpDebug,(int(BRO[0]),int(BRO[1])),1, color = (0,0,0))
                #wpDebug = cv2.circle(wpDebug,(int(BLO[0]),int(BLO[1])),1, color = (0,0,0))
                #wpDebug = cv2.circle(wpDebug,(int(TLI[0]),int(TLI[1])),1, color = (0,0,0))
                #wpDebug = cv2.circle(wpDebug,(int(TRI[0]),int(TRI[1])),1, color = (0,0,0))
                #wpDebug = cv2.circle(wpDebug,(int(BRI[0]),int(BRI[1])),1, color = (0,0,0))
                #wpDebug = cv2.circle(wpDebug,(int(BLI[0]),int(BLI[1])),1, color = (0,0,0))
                #cv2.imshow("WaterPrintOutline", wpDebug)  
                #cv2.waitKey(time)
                #cv2.destroyAllWindows()

        hst = cv2.calcHist([wraped],[0],mask,[256],[0,256])[:int(0.4 * alpha)]
        indexes, _ = scipy.signal.find_peaks(np.array(hst).ravel(), height=20, distance=28)
        if len(indexes) == 1:
            indexes2, _ = scipy.signal.find_peaks(np.array(hst[:80]).ravel(), height=6, distance=10)
            indexes = np.concatenate((indexes, indexes2))
        print(indexes)
        Cut = cutHist(wraped,int(max(indexes) * 0.5 + min(indexes) * 0.5),int(max(indexes) * 2.7), alpha,1,maskOut)
        cv2.imshow("Cut", Cut)  
        cv2.waitKey(time)
        cv2.destroyAllWindows()
        plt.plot(hst,'r')
        plt.xlim([0,int(0.4 * alpha)])
        plt.show()
        wpCnt.extend(waterCnt)

    

    #################### 异常 ################################################################################################################
    if len(screenCnt) < 2 :
        raise Exception("WARNING: IDcard PROCESS FAILED")
    if len(wpCnt) < 1 :
        print("WARNING: WaterPrint PROCESS FAILED: " + imagePath)
        #raise Exception("WARNING: WaterPrint PROCESS FAILED")
    ##=======================================================================================================================================##

    # 展示结果
    if show:
        print("STEP 2: 获取轮廓")
        for sc in screenCnt:
            copy = image.copy()
            cv2.drawContours(copy, [sc], 0, (0, 255, 0), 2)
            cv2.imshow("Outline", copy)
            cv2.waitKey(time)
            cv2.destroyAllWindows()

    #TODO:去水印 : 1.身份证拉正。 2.找水印外围框 3. 水印mask直方图（均衡化） 4.两个峰取阈值 5.原图拉到背景色####################################
        
    #copy = gray.copy()
    #icopy = image.copy()
    #if len(wpCnt) > 0 :
    #    for wc in wpCnt:
    #        waterPrint = four_point_transform(copy, wc)
            #drawHistoGram(waterPrint)
            #Hist = cv2.calcHist([waterPrint],[0],None,[256],[0,256])
            #alpha = 0
            #for i in range(1,256):
            #    m = max(Hist)
            #    if Hist[i] == m:
            #        alpha = i -140
            #        break
            #wpMask = wp * alpha
            #wpWrap = four_point_transform(wpMask, wc, 1, 1000, 1000)
            #backGround = np.ones((1000,1000),np.uint8) * alpha
            #outPut = cv2.add(copy,wpWrap)
            if showprint:
                #print(alpha)
                cv2.imshow("WaterPrint", waterPrint)
                cv2.waitKey(time)
                cv2.destroyAllWindows()
    #wpres = four_point_transform(icopy, waterCnt[1])
    #cv2.imwrite(path + '/waterprint/waterprint.jpg', wpres)


    ############## 透视变换 ###########################################################################################

    
    
    screenCnt = np.array(screenCnt)
    warped1 = four_point_transform(afterProcess, screenCnt[0])
    warped2 = four_point_transform(afterProcess, screenCnt[1])

    # warped = resize(warped, height=600, width=600)
    #warped = cv2.resize(warped, (600, 600))
    if showID:
        cv2.imshow("Scanned1", warped1)
        cv2.waitKey(time)
        cv2.imshow("Scanned2", warped2)
        cv2.waitKey(time)
        cv2.destroyAllWindows()


    #image = resize(orig, height=500)

    return (warped1,warped2)
    #cv2.imwrite('output/transform_ci_zhuan_03.jpg', warped)
