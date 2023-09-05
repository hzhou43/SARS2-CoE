C     Find cros-interface distances in a parallel helix pentamer.
C     NP: number of protomer
C     Note NP can be set to any number and so this code applies to any
C     oligomer.
C     NN: number of residues
C     i,j: residue index; i=1 means residue 11
C     Npairs: maximum number of residue pairs wirh Ca-Ca distances < Dcut
C     KK: actunal number of such pairs
C     k: pair index
C     R: radius of helix (at Ca positions)
C     H: rise per residue;
C     D: spacing b/w helix rims (at Ca posiitons)
C     Dcut: cutoff distance for 13C-13C cross peaks
C     NS: helix rotation in degs, scanned from 10 to 360 deg
C     d2p: degree to radian conversion factor
C     theta1, theta2: helix rotation angles
C
      implicit none

      integer NN,NP,Npairs,NS
      parameter (NN=27,NP=5,Npairs=50)

      real Dcut
      parameter (Dcut=7.2)

      integer VV(5),LL(9),FF(3)
      data VV /14,17,24,25,29/
      data LL /12,18,19,21,27,28,31,34,37/
      data FF /20,23,26/

      real R,H,D
      parameter (R=2.3,H=1.5,D=5.1)

      real x1(NN),y1(NN),z1(NN),x2(NN),y2(NN),z2(NN)
      real pi,d2p,theta0,theta1,theta2
      real dist
      integer II(Npairs),JJ(Npairs)
      real dis(Npairs)

      integer i,j,k,KK,l,m

      pi=acos(-1.0)
      d2p=pi/180.0

      do NS=10,360,10
      theta0=real(NS)*d2p

      theta1=theta0
      theta2=theta0+(360.0/real(NP))*d2p
      do i=1,NN
         x1(i)=R*cos(theta1)
         y1(i)=R*sin(theta1)
         z1(i)=real(i-1)*H

         x2(i)=(2.0*R+D)*cos(180.0/real(NP)*d2p)+R*cos(theta2)
         y2(i)=(2.0*R+D)*sin(180.0/real(NP)*d2p)+R*sin(theta2)

         z2(i)=real(i-1)*H

         theta1=theta1-100.0*d2p
         theta2=theta2-100.0*d2p
      enddo

      k=0
      do i=1,NN
      do j=1,NN
         dist=(x2(i)-x1(j))**2+(y2(i)-y1(j))**2+(z2(i)-z1(j))**2
         dist=dist**0.5
         if (dist .lt. Dcut) then
            k=k+1
            II(k)=i+10
            JJ(k)=j+10
            dis(k)=dist
         endif
      enddo
      enddo
      KK=k

C     filters
C     24-25
      do k=1,KK
         if (II(k) .eq. 24 .and. JJ(k) .eq. 25) goto 500
         if (II(k) .eq. 25 .and. JJ(k) .eq. 24) goto 500
      enddo

C     Phe     
      do k=1,KK
      do l=1,3
         if (II(k) .eq. FF(l) .or. JJ(k) .eq. FF(l)) goto 500
      enddo
      enddo

C     Leu-Val
      do k=1,KK
      do l=1,5
      do m=1,9
         if (II(k) .eq. VV(l) .and. JJ(k) .eq. LL(m)) goto 400
         if (JJ(k) .eq. VV(l) .and. II(k) .eq. LL(m)) goto 400
      enddo
      enddo
      enddo
      goto 500

400   print*,NS
      do k=1,KK
         print*,II(k),JJ(k),dis(k)
      enddo

500   enddo

      end
