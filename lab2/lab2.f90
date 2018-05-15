program lab2
	real, dimension(4) :: array

	array = (/4.0, 3.0, 2.0,1.5 /)

	call make_matrix(5, array)
end program lab2

subroutine make_matrix(n, array_with_values)
	real, dimension(n,n):: matrix
	real, dimension(n-1) :: array_with_values
	
	do i = 1,n
		do j = 1,n
			if (j .le. i) then
				matrix(i,j) = 1
			else
				matrix(i,j) = array_with_values(i)
			end if
		end do
	end do

	do i = 1,n
		do j = 1,n
			write(*,"(f4.1,$)") matrix(i,j)
		end do
		write(*,*)
	end do

end subroutine make_matrix


