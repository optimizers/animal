      program conver
*
*     This program converts the pig breeding matrix data as given
*     originally by Andreas Hofer into the standard Harwell-Boeing Format.
*
*                                  Markus Hegland, 1991
*
      character datfil*6, pedfil*5, daform*22
*     Characterization of input data
*
      parameter(r011=0.0436, r012=0.0141, r022=0.0969)
      parameter(g011=0.0723, g012=0.0105, g022=0.0873)
      parameter(p011=0.0450, p012=0.0393, p022=0.1571)
*
*     First trivial data set:
*
       ! parameter(nme=10, nno=2, ngroup=2)
       ! parameter(nlits=4, nherds=2)
       ! parameter(datfil='x.data',pedfil='x.ped')
       ! parameter(daform='(i5,i2,i4,i5,i4,2f4.0)')
*
*     Small data set:
*
      ! parameter(nme=1181,nno=418,ngroup=6)
      ! parameter(nlits=360,nherds=23)
      ! parameter(datfil='s.data',pedfil='s.ped')
      ! parameter(daform='(i5,i2,i4,i5,i4,2f4.0)')
*
*     Intermediate data set:
*
      ! parameter(nme=3378,nno=1326,ngroup=6)
      ! parameter(nlits=1315,nherds=94)
      ! parameter(datfil='m.data',pedfil='m.ped')
      ! parameter(daform='(i5,i2,i4,i5,i4,2f4.0)')
*
*     Large data set:
*
      ! parameter(nme=11104,nno=1624,ngroup=6)
      ! parameter(nlits=4422,nherds=108)
      ! parameter(datfil='l.data',pedfil='l.ped')
      ! parameter(daform='(i5,i2,i4,i5,i4,2f4.0)')
*
*     Very large data set:
*
      parameter(nme=70550,nno=6295,ngroup=6)
      parameter(nlits=26798,nherds=2233)
      parameter(datfil='v.data',pedfil='v.ped')
      parameter(daform='(i5,i2,i4,i5,i4,2f4.0)')
*
*     Extreme data set: (different daform)
*
      ! parameter(nme=197949,nno=14188,ngroup=6)
      ! parameter(nlits=74240,nherds=4964)
      ! parameter(datfil='e.data',pedfil='e.ped')
      ! parameter(daform='(i6,i2,i4,i5,i4,2f4.0)')
*
      parameter(n=2*(nme+nno+ngroup+nlits+nherds))
      parameter(m=2*(nme+nlits+nme+nno))
      parameter(nm=3*(nme*3+(nme+nno)*3+nlits))
*
      integer rowptr(m+1), colind(nm), colptr(n+1), rowind(nm)
      integer ani, sex, herd, litter, t, dad, mum, type
      double precision values(nm), rhs(m), valcol(nm)

      character title*72, key*8, mxtype*3,
     *          ptrfmt*16, indfmt*16, valfmt*20, rhsfmt*20

      data title/'PIG BREEDING MATRIX 2. A.HOFER/M.HEGLAND, 1991'/
      data key/'yyy'/
      data mxtype/'RRA'/
      data ptrfmt/'(12I8)'/, indfmt/'(12I8)'/  
      data valfmt/'(4D22.16)'/, rhsfmt/'(4D22.16)'/
*
      data rhs/m*0.D0/
*
*     Input data
*
      open(1,file=datfil,status='old')
      ic = 1
      vc11 = r011
      vc12 = r012
      vc22 = r022
      do 1 i=1,2*nme,2
        rowptr(i) = ic
        rowptr(i+1) = ic + 6
        read(1,daform) ani,sex,herd,litter,t,y,yprime
        rhs(i) = y
        rhs(i+1) = yprime
        colind(ic)  = 2*ani-1
        colind(ic+1) = 2*ani
        colind(ic+6) = 2*ani
        colind(ic+2) = 2*(nme+nno+litter)-1
        colind(ic+3) = 2*(nme+nno+litter)
        colind(ic+7) = 2*(nme+nno+litter)
        colind(ic+4) = 2*(nme+nno+nlits+herd)-1
        colind(ic+5) = 2*(nme+nno+nlits+herd)
        colind(ic+8) = 2*(nme+nno+nlits+herd)
        values(ic) = vc11
        values(ic+1) = vc12
        values(ic+6) = vc22
        values(ic+2) = vc11
        values(ic+3) = vc12
        values(ic+7) = vc22
        values(ic+4) = vc11
        values(ic+5) = vc12
        values(ic+8) = vc22
        ic = ic+9
1     continue
*
      open(2,file=pedfil,status='old')
      do 2 i=2*nme+1,2*(2*nme+nno),2
        read(2,'(3i6,i2,f13.10)') ani, dad, mum, type, fac
        if (dad.ne.mum) then
          ioff = 6
        else
          ioff = 4
        endif
        rowptr(i) = ic
        rowptr(i+1) = ic+ioff
        rhs(i) = 0.D0
        rhs(i+1) = 0.D0
        vc11 = g011/sqrt(fac)
        vc12 = g012/sqrt(fac)
        vc22 = g022/sqrt(fac)
        colind(ic) = 2*ani-1
        colind(ic+1) = 2*ani
        colind(ic+ioff) = 2*ani
        values(ic) = vc11
        values(ic+1) = vc12
        values(ic+ioff) = vc22
        if (type.eq.3) then
          colind(ic+2) = 2*dad-1
          colind(ic+3) = 2*dad
          colind(ic+ioff+1) = 2*dad
          colind(ic+4) = 2*mum-1
          colind(ic+5) = 2*mum
          colind(ic+ioff+2) = 2*mum
          values(ic+2) = -vc11/2.0
          values(ic+3) = -vc12/2.0
          values(ic+ioff+1) = -vc22/2.0
          values(ic+4) = -vc11/2.0
          values(ic+5) = -vc12/2.0
          values(ic+ioff+2) = -vc22/2.0
          ic = ic+ioff+3
        elseif (type.eq.1) then
          if (dad.ne.mum) then
            colind(ic+2) = 2*(dad + nlits + nherds)-1
            colind(ic+3) = 2*(dad + nlits + nherds)
            colind(ic+ioff+1) = 2*(dad + nlits + nherds)
            colind(ic+4) = 2*(mum + nlits + nherds)-1
            colind(ic+5) = 2*(mum + nlits + nherds)
            colind(ic+ioff+2) = 2*(mum + nlits + nherds)
            values(ic+2) = -vc11/2.0
            values(ic+3) = -vc12/2.0
            values(ic+ioff+1) = -vc22/2.0
            values(ic+4) = -vc11/2.0
            values(ic+5) = -vc12/2.0
            values(ic+ioff+2) = -vc22/2.0
            ic = ic+ioff+3
          else
            colind(ic+2) = 2*(dad + nlits + nherds)-1
            colind(ic+3) = 2*(dad + nlits + nherds)
            colind(ic+ioff+1) = 2*(dad + nlits + nherds)
            values(ic+2) = -vc11
            values(ic+3) = -vc12
            values(ic+ioff+1) = -vc22
            ic = ic+ioff+2
          endif 
        else
          if (dad.gt.nme+nno) then
            colind(ic+2) = 2*(dad + nlits + nherds)-1
            colind(ic+3) = 2*(dad + nlits + nherds)
            colind(ic+ioff+1) = 2*(dad + nlits + nherds)
          else
            colind(ic+2) = 2*dad-1
            colind(ic+3) = 2*dad
            colind(ic+ioff+1) = 2*dad
          endif
          if (mum.gt.nme+nno) then
            colind(ic+4) = 2*(mum + nlits + nherds)-1
            colind(ic+5) = 2*(mum + nlits + nherds)
            colind(ic+ioff+2) = 2*(mum + nlits + nherds)
          else
            colind(ic+4) = 2*mum-1
            colind(ic+5) = 2*mum
            colind(ic+ioff+2) = 2*mum
          endif
          values(ic+2) = -vc11/2.0
          values(ic+3) = -vc12/2.0
          values(ic+ioff+1) = -vc22/2.0
          values(ic+4) = -vc11/2.0
          values(ic+5) = -vc12/2.0
          values(ic+ioff+2) = -vc22/2.0
          ic = ic+ioff+3
        endif
2     continue
*
      vc11 = p011
      vc12 = p012
      vc22 = p022
      do 3 i=2*(2*nme+nno)+1,2*(2*nme+nno+nlits),2
        rowptr(i) = ic
        rowptr(i+1) = ic+2
        colind(ic) = i-2*nme-1
        colind(ic+1) = i-2*nme
        colind(ic+2) = i-2*nme
        rhs(i) = 0.D0
        rhs(i+1) = 0.D0
        values(ic) = vc11
        values(ic+1) = vc12
        values(ic+2) = vc22
        ic = ic+3
3     continue
      rowptr(m+1)=ic
      ic=ic-1
*
*     construct the matrix by column (colptr, colind)
*
*Set to zero the column pointers of A.

          do 220 j=1,n+1
             colptr(j)=0
 220      continue

* Count the elements in each column.

          dO 225 k=1,ic
              j=colind(k)
             colptr(j) = colptr(j) +1
 225      continue

* Prepare the column pointers in such a way that pointer colptr(j)
* points to the first position of column J+1.
         colptr(1)=colptr(1)+1
         do 227 k=2,n+1
             colptr(K)= colptr(k)+colptr(k-1)
 227     continue
* Set the column pointers and put the matrix in column order.
         do 230 i=1,m
           i1=rowptr(i)
           i2=rowptr(i+1)-1
           do 235 k=i1,i2
              j=colind(k)
            jpos=colptr(j)
            rowind(jpos-1)=i
            valcol(jpos-1)=values(k)
            colptr(j)=jpos-1
 235       continue
 230    continue

*
*     output data in converted form
*
      npt=(m+11)/12
      nro=(ic+11)/12
      nva=(ic+3)/4
      nrh=(m+3)/4
*
      write(*,'(A72,A8)') title, key
      write(*,'(5I14)') npt+nro+nva+nrh,npt,nro,nva,nrh
      write(*,'(A3,11X,4I14)') mxtype,m,n,ic,0
      write(*,'(2A16,2A20)') ptrfmt, indfmt, valfmt, rhsfmt
      write(*,'(A1,13X,2I14)') 'F',1,m
*
      write(*,ptrfmt) colptr
      write(*,indfmt) (rowind(i),i=1,ic)
      write(*,valfmt) (valcol(i),i=1,ic)
      write(*,rhsfmt) rhs
*
      end
