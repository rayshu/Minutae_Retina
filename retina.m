close all;
im=imread('Retina.tif');
bw_mask=imread('mask.gif');
bw_mask=logical(bw_mask);
figure, imshow(im), title('Green channel retinal image');
im=im(:,:,2);
im=mat2gray(im).*mat2gray(bw_mask);
im=imcomplement(im);
im=im2double(im);

DEG_NUM=12;
LEN_c=11;
LEN_o=11;
LEN_diff=7;

ic1=reconstruction_by_dilation(im,LEN_c,DEG_NUM);
io1=min_openings(im,LEN_o,DEG_NUM);
iv=mat2gray(ic1-io1);
imDiff=smooth_cross_section(iv,LEN_diff,LEN_c);
imL=reconstruction_by_dilation(imDiff,LEN_c,DEG_NUM);
imF=reconstruction_by_erosion(imL,LEN_c,DEG_NUM);
figure,imshow(iv);title('Enhanced blood vessels');
figure,imshow(imL);title('Dilation');
figure,imshow(imF);title('Erosion');
figure,imshow(imDiff);title('Difference');

TH_LOW=30;
TH_HIGH=40;
min_obj=180;
min_hole=10;

mask=im2bw(imF,TH_LOW/255);
marker=im2bw(imF,TH_HIGH/255);
bw_result=imreconstruct(marker,mask);

bw_result=bw_result& bw_mask;
bw_result = clear_bw(bw_result, min_obj, min_hole);
figure,imshow(bw_result);title('Segmented blood vessel');

I=imread('thinn.jpg');
J=I(:,:,1)>100;
K=bwmorph(J,'thin','inf');
figure, imshow(~K), title('Thinned blood vessel');

fun=@minutie;
L = nlfilter(K,[3 3],fun);

%%%%%%%%Termination red color
LTerm=(L==1);
imshow(LTerm)
LTermLab=bwlabel(LTerm);
propTerm=regionprops(LTermLab,'Centroid');
CentroidTerm=round(cat(1,propTerm(:).Centroid));
figure, imshow(~K), title('Minutae considering spurs');
hold on
plot(CentroidTerm(:,1),CentroidTerm(:,2),'ro')

%%%%%%Bifurcation green color
LBif=(L==3);
LBifLab=bwlabel(LBif);
propBif=regionprops(LBifLab,'Centroid','Image');
CentroidBif=round(cat(1,propBif(:).Centroid));
plot(CentroidBif(:,1),CentroidBif(:,2),'go')

%%%remove spurs
D=5;
Distance=DistEuclidian(CentroidBif,CentroidTerm);
SpuriousMinutae=Distance<D;
[i,j]=find(SpuriousMinutae);
CentroidBif(i,:)=[];
CentroidTerm(j,:)=[];

Distance=DistEuclidian(CentroidBif);
SpuriousMinutae=Distance<D;
[i,j]=find(SpuriousMinutae);
CentroidBif(i,:)=[];

Distance=DistEuclidian(CentroidTerm);
SpuriousMinutae=Distance<D;
[i,j]=find(SpuriousMinutae);
CentroidTerm(i,:)=[];
hold off
imshow(~K)
hold on
plot(CentroidTerm(:,1),CentroidTerm(:,2),'ro')
plot(CentroidBif(:,1),CentroidBif(:,2),'go')
hold off

Kopen=imclose(K,strel('square',7));

KopenClean= imfill(Kopen,'holes');
KopenClean=bwareaopen(KopenClean,5);
imshow(KopenClean)
KopenClean([1 end],:)=0;
KopenClean(:,[1 end])=0;
ROI=imerode(KopenClean,strel('disk',10));
figure, imshow(ROI), title('ROI');

imshow(I)
hold on
figure, imshow(K), title('Minutiae after ignoring spurs');

hold on
plot(CentroidTerm(:,1),CentroidTerm(:,2),'ro')
plot(CentroidBif(:,1),CentroidBif(:,2),'go')
hold off

[m,n]=size(I(:,:,1));
indTerm=sub2ind([m,n],CentroidTerm(:,1),CentroidTerm(:,2));
Z=zeros(m,n);
Z(indTerm)=1;
ZTerm=Z.*ROI';
[CentroidTermX,CentroidTermY]=find(ZTerm);

indBif=sub2ind([m,n],CentroidBif(:,1),CentroidBif(:,2));
Z=zeros(m,n);
Z(indBif)=1;
ZBif=Z.*ROI';
[CentroidBifX,CentroidBifY]=find(ZBif);




hold on
plot(CentroidTermX,CentroidTermY,'ro','linewidth',2)
plot(CentroidBifX,CentroidBifY,'go','linewidth',2)