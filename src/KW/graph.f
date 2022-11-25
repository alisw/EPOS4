C
C  This file is part of EPOS4
C  Copyright (C) 2022 research institutions and authors (See CREDITS file)
C  This file is distributed under the terms of the GNU General Public License version 3 or later
C  (See COPYING file for the text of the licence)
C

c#######################################################################
      subroutine RootCanvas(npad)
c#######################################################################
      !   call RootCanvas(npad)
      !   call RootPave(text)
      !   call RootPadMult(npad)
      !   do i=1,npad
      !     call RootPad(i)
      !     call RootHisto(ndim,text,nx,xmin,xmax,ny,ymin,ymax)
      !     do ...
      !       call RootFill(ndim,x,y,z)
      !     enddo
      !     call RootDraw(style,textx,texty)
      !   enddo
      !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#include "aaa.h"
      character cc*5, txt*80
      data icanvas /0/
      save icanvas
      icanvas=icanvas+1
      ixoff=mod(20*(icanvas-1),900)
      iyoff=mod(20*(icanvas-1),450)
      if(icanvas.gt.1)iyoff=iyoff+10
      if(npad.eq.-2)then
      i1=10+ixoff
      i2=10+iyoff
      i3=500
      i4=500
      elseif(npad.eq.1)then
      i1=10+ixoff
      i2=10+iyoff
      i3=700
      i4=500
      elseif(npad.eq.-1)then
      i1=10+ixoff
      i2=10+iyoff
      i3=900
      i4=500
      elseif(npad.eq.-3)then
      i1=10+ixoff
      i2=10+iyoff
      i3=700
      i4=700
      iTheta=60
      iPhi=70
      elseif(npad.eq.-4)then
      i1=10+ixoff
      i2=10+iyoff
      i3=370
      i4=370
      elseif(npad.eq.6)then
      i1=10+ixoff
      i2=10+iyoff
      i3=900
      i4=700
      elseif(npad.eq.8)then
      i1=0
      i2=0
      i3=1200
      i4=700
      endif
      ii=0
      if(icanvas.eq.1)then
      nopenr=nopenr+1
      ii=index(fnhi(1:nfnhi),".histo")-1
      open(95,file=fnhi(1:ii)//".cc")
      write(95,'(12(a/),a,a,4(a/),a)')
      !-----------------------------------------------------------------
     .'{                                                          ',
     .'  gROOT->Reset();                                          ',
     .'  gROOT->SetStyle("Plain");                                ',
     .'  gStyle->SetPalette(1);                                   ',
     .'  gStyle->SetCanvasColor(10);   //33=grey 18=light grey    ',
     .'  gStyle->SetFrameFillColor(10);   //10=white              ',
     .'  const Int_t NRGBs = 6;                                   ',
     .'  const Int_t NCont = 20;                                  ',
     .'  Double_t stops[NRGBs] = {0.0,0.05,0.30,0.70,0.85,1.00 }; ',
     .'  Double_t red[NRGBs]   = {0.0,0.00,0.00,1.00,1.00,1.00 }; ',
     .'  Double_t green[NRGBs] = {0.0,0.00,0.81,1.00,0.50,0.00 }; ',
     .'  Double_t blue[NRGBs]  = {0.5,1.00,0.80,0.00,0.00,0.00 }; ',
     .'  TColor::CreateGradientColorTable',
     .                   '(NRGBs, stops, red, green, blue, NCont);',
     .'  gStyle->SetNumberContours(NCont);                   ',
     .'  gStyle->SetPadBottomMargin(0.1);                    ',
     .'  gStyle->SetPadRightMargin(0.12);                    ',
     .'  gStyle->SetPadLeftMargin(0.1);                      ',
     .'  gStyle->SetPadTopMargin(0.1);                       '
      !-----------------------------------------------------------------
c     .'  Double_t stops[NRGBs] = {0,0.15,0.50,0.80,0.90,1.00 }; ',
c     .'  Double_t red[NRGBs]   = {1,0.00,0.00,1.00,1.00,1.00 }; ',
c     .'  Double_t green[NRGBs] = {1,0.00,0.81,1.00,0.80,0.00 }; ',
c     .'  Double_t blue[NRGBs]  = {1,1.00,0.80,0.00,0.00,0.00 }; ',
      !-----------------------------------------------------------------
      endif
      write(cc,'(i5)')icanvas
      n=5-(int(log10(float(icanvas)))+1)
      do m=2,n
      cc(m:m)='0'
      enddo
      cc(1:1)='_'
      ii1=0
      do i=1,ii
      if(fnhi(i:i).eq.'/') ii1=i
      enddo
      ii1=ii1+1
      iii=ii-ii1+1
      txt(1:iii)=fnhi(ii1:ii)
      do m=1,iii
      if(txt(m:m).eq.'-')txt(m:m)='_'
      enddo
      write(95,'(a,a,4(i4,a)/a)')
      !----------------------------------------------------------------'
     .'  TCanvas *'//txt(1:iii)//cc//' = new TCanvas("',
     .   txt(1:iii)//cc//'", "'//txt(1:iii)//cc//'",',
     .   i1,',',i2,',',i3,',',i4,');',
     .'  '//txt(1:iii)//cc//'->Range(0,0,25,18);           '
      !----------------------------------------------------------------'
      if(npad.eq.-3)write(95,'(a,i3,a)')
      !----------------------------------------------------------------'
     . '  '//txt(1:iii)//cc//'->SetTheta(',iTheta,');           '
     .,'  '//txt(1:iii)//cc//'->SetPhi(',iPhi,');           '
      !----------------------------------------------------------------'
      end

c#######################################################################
      subroutine RootPave(text)
c#######################################################################
      character text*(*), cc*5
      data ipave /0/
      save ipave
      ipave=ipave+1
      write(cc,'(i5)')ipave
      n=4-ipave/10
      cc(1:1)='p'
      do m=2,n
      cc(m:m)='0'
      enddo
      kmax=1
      do while(text(kmax:kmax).ne.';')
      kmax=kmax+1
      enddo
      kmax=kmax-1
      write(95,'(2a/4(a/))')
      !-----------------------------------------------------------------
     .'  TPaveLabel *'//cc//' = new TPaveLabel(1,16.3,24,17.5,"'    ,
     .         text(1:kmax)//'","br");'                             ,
     .'  '//cc//'->SetFillColor(18);                               ',
     .'  '//cc//'->SetTextFont(32);                                ',
     .'  '//cc//'->SetTextColor(49);                               ',
     .'  '//cc//'->Draw();                                         '
      !-----------------------------------------------------------------
      end

c#######################################################################
      subroutine RootPadMult(n)
c#######################################################################
      if(n.eq.8)then
      write(95,'(16(a/))')
      !-----------------------------------------------------------------
     .'  pad1 = new TPad("pad1","",0.01,0.45,0.24,0.89,33);  ',
     .'  pad2 = new TPad("pad2","",0.25,0.45,0.49,0.89,33);  ',
     .'  pad3 = new TPad("pad3","",0.50,0.45,0.74,0.89,33);  ',
     .'  pad4 = new TPad("pad4","",0.75,0.45,0.99,0.89,33);  ',
     .'  pad5 = new TPad("pad5","",0.01,0.01,0.24,0.44,33);  ',
     .'  pad6 = new TPad("pad6","",0.25,0.01,0.49,0.44,33);  ',
     .'  pad7 = new TPad("pad7","",0.50,0.01,0.74,0.44,33);  ',
     .'  pad8 = new TPad("pad8","",0.75,0.01,0.99,0.44,33);  ',
     .'  pad1->Draw();                                       ',
     .'  pad2->Draw();                                       ',
     .'  pad3->Draw();                                       ',
     .'  pad4->Draw();                                       ',
     .'  pad5->Draw();                                       ',
     .'  pad6->Draw();                                       ',
     .'  pad7->Draw();                                       ',
     .'  pad8->Draw();                                       '
      !-----------------------------------------------------------------
      elseif(n.eq.6)then
      write(95,'(12(a/))')
      !-----------------------------------------------------------------
     .'  pad1 = new TPad("pad1","",0.010,0.45,0.323,0.89,33);  ',
     .'  pad2 = new TPad("pad2","",0.333,0.45,0.656,0.89,33);  ',
     .'  pad3 = new TPad("pad3","",0.666,0.45,0.989,0.89,33);  ',
     .'  pad4 = new TPad("pad4","",0.010,0.01,0.323,0.44,33);  ',
     .'  pad5 = new TPad("pad5","",0.333,0.01,0.656,0.44,33);  ',
     .'  pad6 = new TPad("pad6","",0.666,0.01,0.989,0.44,33);  ',
     .'  pad1->Draw();                                       ',
     .'  pad2->Draw();                                       ',
     .'  pad3->Draw();                                       ',
     .'  pad4->Draw();                                       ',
     .'  pad5->Draw();                                       ',
     .'  pad6->Draw();                                       '
      !-----------------------------------------------------------------
      endif
      end

c#######################################################################
      subroutine RootPad(i)
c#######################################################################
      write(95,'(2(a,i1,a/))')
      !-----------------------------------------------------------------
     .'  pad',i,'->cd();                        ',
     .'  pad',i,'->GetFrame()->SetFillColor(19);'
      !-----------------------------------------------------------------
      end

c#######################################################################
      subroutine RootHisto(ndim,text,nx,xmin,xmax,ny,ymin,ymax)
c#######################################################################
      character text*(*)
      character xx*5,yy*5,tt*5,ttit*80
      common /canvas/xx,yy,tt,ttit
      data ihisto /0/
      save ihisto
      ihisto=ihisto+1
      ttit=text
      write(xx,'(i5)')ihisto
      xx(1:1)='x'
      n=4-int(log10(float(ihisto)))
      do m=2,n
      xx(m:m)='0'
      enddo
      yy=xx
      yy(1:1)='y'
      tt=xx
      tt(1:1)='t'

      if(ndim.eq.2)
     . write(95,'(2a,2(i5,a,e11.3,a,e11.3,a)/a,i5,a/a,i5,a)')
      !-----------------------------------------------------------------
     .   '  '//xx//' = new TH2F("'//xx//'","','",'
     .          ,nx,',',xmin,',',xmax,',',ny,',',ymin,',',ymax,');'
     . , '  '//xx//'->GetXaxis()->SetRange(1,',nx,');'
     . , '  '//xx//'->GetYaxis()->SetRange(1,',ny,');'
      !-----------------------------------------------------------------
      if(ndim.eq.1)
     . write(95,'(2a,2(i5,a,e11.3,a,e11.3,a)/a,i5,a)')
      !-----------------------------------------------------------------
     .   '  '//xx//' = new TH1F("'//xx//'","','",'
     .          ,nx,',',xmin,',',xmax,');'
     . , '    '//xx//'->GetXaxis()->SetRange(1,',nx,');'
      !-----------------------------------------------------------------
      end

c#######################################################################
      subroutine RootFill(ndim,x,y,z)
c#######################################################################
      character xx*5,yy*5,tt*5,ttit*80
      common /canvas/xx,yy,tt,ttit
      if(ndim.eq.2)write(95,'(a,3(e11.3,a))')
      !-----------------------------------------------------------------
     . '  '//xx//'->Fill(',x,',',y,',',z,');'
      !-----------------------------------------------------------------
      if(ndim.eq.1)write(95,'(a,2(e11.3,a))')
      !-----------------------------------------------------------------
     . '  '//xx//'->Fill(',x,',',y,');'
      !-----------------------------------------------------------------
      end

c#######################################################################
      subroutine RootDraw(style,textx,texty)
c#######################################################################
      ! draw styles:  lego1...4 surf1...4 contz colz contz&colz
      !........................................................
      character text*80,style*(*),textx*(*),texty*(*),cinf*45
      character stylel*256,textxl*256,textyl*256
      character coffx*4, coffy*4
      character xx*5,yy*5,tt*5,ttit*80
      common /canvas/xx,yy,tt,ttit
      text=ttit
      lmax=1
      do while(style(lmax:lmax).ne.';')
      lmax=lmax+1
      enddo
      lmax=min(256,lmax-1)
      stylel=style(1:lmax)
      kx=1
      do while(textx(kx:kx).ne.';')
      kx=kx+1
      enddo
      kx=min(256,kx-1)
      textxl=textx(1:kx)
      ky=1
      do while(texty(ky:ky).ne.';')
      ky=ky+1
      enddo
      ky=min(256,ky-1)
      textyl=texty(1:ky)
      cinf='//  lego1...4 surf1...4 contz colz contz&colz'
      if(style(1:lmax).eq.'cont4Z')then
        coffx='1.00'
        coffy='1.00'
      elseif(style(1:lmax).eq.'contz')then
        coffx='1.00'
        coffy='1.00'
      elseif(style(1:lmax).eq.'colz')then
        coffx='1.00'
        coffy='1.15'
      else
        coffx='2.00'
        coffy='2.00'
      endif
      kmax=1
      do while(text(kmax:kmax).ne.';')
      kmax=kmax+1
      enddo
      kmax=kmax-1

      write(95,'(a,a/15(a/))')
      !-----------------------------------------------------------------
     .'  '//xx//'->SetMinimum('//xx//'->GetMinimum()*1.05-'
     .           ,xx//'->GetMaximum()*0.05);',
     .'  '//xx//'->SetStats(0);                    ',
     .'  '//xx//'->Draw("'//stylel(1:lmax)//'");    '//cinf,
     .'  '//xx//'->SetMarkerStyle(20) ;            ',
     .'  '//xx//'->SetMarkerColor(2) ;             ',
     .'  '//xx//'->SetMarkerSize(0.4) ;            ',
     .'  '//xx//'->GetXaxis()->SetTitleOffset('//coffx//');',
     .'  '//xx//'->SetXTitle("'//textxl(1:kx)//'"); ',
     .'  '//xx//'->GetYaxis()->SetTitleOffset('//coffy//');',
     .'  '//xx//'->SetYTitle("'//textyl(1:ky)//'"); ',
     .'  '//yy//' = new TPaveText(0.5,0.88,0.5,0.97,"NDC");',
     .'  '//yy//'->SetFillColor(0);',
     .'  '//yy//'->SetBorderSize(0);',
     .'  '//tt//' = '//yy//'->AddText(1.,1.,"'//text(1:kmax)//'");',
     .'  '//tt//'->SetTextSize(0.036);',
     .'  '//yy//'->Draw();'
      !-----------------------------------------------------------------
      end

c#######################################################################
      subroutine RootClose()
c#######################################################################
#include "aaa.h"
      if(nopenr.gt.0)then
      ii=index(fnhi(1:nfnhi),".histo")-1
      write(95,'(a)')'}'
      close (95)
      endif
      end


