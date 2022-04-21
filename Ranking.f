      module Spop_Ranking
********************************************************************************
c     ranks indivuduals.
********************************************************************************
      use GAparam, only: mobj,mcon
      implicit none

*   work variables, frequently used in this module:
      integer:: iobj,icn,ipop,isub

      private
      public Ranking
      contains



      subroutine Ranking (iobj,icon,n,f,g,v,a,b,mrank,nrank,nrk)
********************************************************************************
c     ranks individuals.
c
c     subprograms used:
c       Violate ......... calculates constraint violations.
c       ParetoRanking ... performs Pareto ranking.
c       FrontRanking .... performs Pareto front ranking.
c
c     module used:
        use GAparam, only: mprank
********************************************************************************
      real:: f(mobj,n),g(mcon,n),v(mcon,n),a(mobj,mobj),b(mobj)
      integer:: n,mrank(n),nrank(n),nrk(n)
      integer:: iobj(2),icon(2)

*  setup:
      call Violate(icon(1),icon(2),n,g,v)

*   ranking:
      if (mprank == 1) call ParetoRanking(iobj(1),iobj(2),icon(1),icon(2),n,f,g,v,a,b,mrank,nrank,nrk)
      if (mprank == 2) call  FrontRanking(iobj(1),iobj(2),icon(1),icon(2),n,f,g,v,a,b,mrank,nrank,nrk)

      return
      end subroutine Ranking



      subroutine Violate (icon1,icon2,n,g,v)
********************************************************************************
c     calculates violations to constraints.
********************************************************************************
      integer:: n, ip, icon, icon1, icon2
      real:: g(mcon,n), v(mcon,n)

      do ip = 1, n
        do icon = icon1, icon2
          v(icon,ip)= max(0.,g(icon,ip))
        enddo
      enddo

      return
      end subroutine Violate



      subroutine ParetoRanking (iobj1,iobj2,icon1,icon2,n,f,g,v,a,b,mrank,nrank,nrk)
********************************************************************************
c	  ranks indivisuals by pareto ranking.
********************************************************************************
      real:: f(mobj,n),g(mcon,n),v(mcon,n),a(mobj,mobj),b(mobj)
      integer:: mrank(n),nrank(n),nrk(n)
      integer:: n,iobj1,iobj2,icon1,icon2,irank

*   setup:
      mrank(:)= 0; nrk(:)= 1

*   ranking:
      do ipop = 1, n
        do isub = 1, n
          call RankCounter(iobj1,iobj2,icon1,icon2,f(1,ipop),f(1,isub),v(1,ipop),v(1,isub),a,b,nrk(ipop))
        enddo
      enddo

*   rank, pareto:
      do ipop = 1, n
        irank = nrk(ipop)
        nrank(ipop) = irank
        mrank(irank) = mrank(irank)+1
      enddo

      return
      end subroutine ParetoRanking



      subroutine FrontRanking (iobj1,iobj2,icon1,icon2,n,f,g,v,a,b,mrank,nrank,nrk)
********************************************************************************
c     ranks indivisuals by pareto front ranking.
********************************************************************************
      real:: f(mobj,n),g(mcon,n),v(mcon,n),a(mobj,mobj),b(mobj)
      integer:: mrank(n),nrank(n),nrk(n)
      integer:: n,iobj1,iobj2,icon1,icon2,irank,inum

*   setup:
      mrank(:)= 0; nrk(:)= 1; inum = 1; irank = 0

      do while (inum <= n)

*   ranking:
        do ipop = 1, n
          if (nrk(ipop) < 0) cycle
          nrk(ipop) = 1
          do isub = 1, n
            if (nrk(isub) < 0) cycle

            call RankCounter(iobj1,iobj2,icon1,icon2,f(1,ipop),f(1,isub),v(1,ipop),v(1,isub),a,b,nrk(ipop))
          enddo
        enddo

*   rank, pareto front:
        irank = irank + 1
        do ipop = 1, n
          if (nrk(ipop) == 1) then
            nrk(ipop) = -1
            inum = inum + 1
            nrank(ipop) = irank
            mrank(irank) = mrank(irank) + 1
          endif
        enddo

      enddo

      return
      end subroutine FrontRanking



      subroutine RankCounter (iobj1,iobj2,icon1,icon2,fp,fs,vp,vs,a,b,nrk)
********************************************************************************
c     counts up rank of an individual with constraints and objectives.
********************************************************************************
      integer:: iobj1,iobj2,icon1,icon2,nrk,np,ns,npx,nsx,idomC,idomO
      real:: fp(mobj),fs(mobj),a(mobj,mobj),b(mobj)
      real:: vp(mcon),vs(mcon)

*   More Less-Violations Method:
      call CompareC(icon1,icon2,vp,vs,np,ns,npx,nsx,idomC)

      if (npx > 0 .or. nsx > 0) then
         if (idomC == 1) nrk=nrk+1
      else
         call CompareO(iobj1,iobj2,fp,fs,a,b,idomO)

         if (np > 0 .or. ns > 0) then
            if (idomC == 1 .and. idomO == 1) nrk=nrk+1
         else
            if (idomO == 1) nrk=nrk+1
         endif
      endif

      return
      end subroutine RankCounter



      subroutine CompareC (icon1,icon2,vp,vs,npop,nsub,npopx,nsubx,idom)
********************************************************************************
c     examines whether or not individual 'pop' is dominated by individual 'sub' with respect to constraints.
c
c     idom ... flag that indicates domination.
c              =-1 not dominated
c              = 0 identical
c              = 1 dominated
c
c     module used:
               use GAparam, only: gw, vx
********************************************************************************
      real:: vp(mcon), vs(mcon), vpop, vsub, tshd, Fpop, Fsub
      integer:: icon1, icon2, icon, idom, npop, nsub, npopx, nsubx

         Fpop = 0.
         Fsub = 0.
         npop = 0
         nsub = 0
         npopx = 0
         nsubx = 0

         do icon = icon1, icon2
            vpop = vp(icon)
            vsub = vs(icon)
            tshd = vx(icon)

            if ( vpop > tshd ) npopx = npopx + 1
            if ( vsub > tshd ) nsubx = nsubx + 1
            if ( vpop > 0.   ) npop = npop + 1
            if ( vsub > 0.   ) nsub = nsub + 1
            if ( vpop > vsub ) Fpop = Fpop + gw(icon)
            if ( vsub > vpop ) Fsub = Fsub + gw(icon)
         enddo

              if (Fpop <  Fsub) then; idom =-1
         else if (Fpop == Fsub) then; idom = 0
         else                       ; idom = 1 ; endif

      return
      end subroutine CompareC



      subroutine CompareO (iobj1,iobj2,fpop,fsub,a,b,idom)
********************************************************************************
c     examines whether individual 'pop' is dominated by individual 'sub' with respect to objectives.
c
c     idom ... flag that indicates that the domination exists.
c              =-1 not dominated
c              = 0 identical
c              = 1 dominated
c
********************************************************************************
      implicit none
      real:: fpop(mobj),fsub(mobj),a(mobj,mobj),b(mobj),df(mobj),sum,d
      integer:: iobj1,iobj2,idom,i,j

*   difference, normalized:
      do j = iobj1, iobj2
         df(j) = (fpop(j) - fsub(j))/ b(j)
      enddo

*   dominance check:
      idom = 1
      sum  = 0.
      do i = iobj1, iobj2
         d = 0.
         do j = iobj1, iobj2
            d = d + a(i,j)* df(j)
         enddo
         if (d > 0.) then
            idom = -1
            return
         else
            sum = sum + d
         endif
      enddo
      if (sum == 0.) idom = 0
      return
      end subroutine CompareO



      end module Spop_Ranking
