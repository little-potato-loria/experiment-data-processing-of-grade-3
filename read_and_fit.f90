program read_and_fit
    implicit none
    integer, parameter :: max_rows = 10000
    integer :: i, num_rows, ios
    real(8), allocatable :: x_data(:), y_data(:)
    real(8) :: A, B, C, chisq_old, chisq_new
    real(8), dimension(3) :: delta
    character(len=256) :: line
    character(len=256) :: filename

    ! 从命令行获取文件名
    call get_command_argument(1, filename)
    if (len_trim(filename) == 0) then
        print *, 'Error: No filename provided.'
        stop
    end if

    ! 预先分配内存给数组
    allocate(x_data(max_rows), stat=ios)
    if (ios /= 0) then
        print *, 'Error: Unable to allocate memory for x_data.'
        stop
    end if
    allocate(y_data(max_rows), stat=ios)
    if (ios /= 0) then
        print *, 'Error: Unable to allocate memory for y_data.'
        stop
    end if

    ! 打开文件
    open(unit=10, file=filename, status='old', action='read', iostat=ios)
    if (ios /= 0) then
        print *, 'Error: Unable to open file ', trim(filename)
        stop
    end if

    ! 初始化行数
    num_rows = 0

    ! 读取文件中的每一行
    do
        read(10, '(A)', iostat=ios) line
        if (ios /= 0) exit
        num_rows = num_rows + 1
        if (num_rows > max_rows) then
            print *, 'Error: Too many rows in the CSV file.'
            exit
        end if

       ! 从行中提取数据
        read(line, *) x_data(num_rows), y_data(num_rows)
        y_data(num_rows) = y_data(num_rows)**2  ! 将第二列数据平方
    end do

    ! 关闭文件
    close(10)

    ! 输出前3行数据
    print *, 'First 3 rows of data:'
    do i = 1, min(num_rows, 3)
        print *, 'Row ', i, ': x = ', x_data(i), ', y = ', y_data(i)
    end do

    ! 初始化拟合参数
    A = 5E-6
    B = -0.007
    C = 3
    chisq_old = 1.0e10

    ! 执行拟合
    call lmfit(x_data, y_data, num_rows, A, B, C, chisq_old)

    ! 输出拟合参数
    print *, 'Fitted parameters:'
    print *, 'A = ', A
    print *, 'B = ', B
    print *, 'C = ', C

contains

    function quadratic(x, A, B, C) result(y)
        real(8), intent(in) :: x, A, B, C
        real(8) :: y
        y = A * x**2 + B * x + C
    end function quadratic

    subroutine lmfit(x, y, n, A, B, C, chisq_old)
        implicit none
        integer, intent(in) :: n
        real(8), dimension(n), intent(in) :: x, y
        real(8), intent(inout) :: A, B, C, chisq_old
        real(8), dimension(n, 3) :: J
        real(8), dimension(n) :: L, res
        real(8), dimension(3, 3) :: H
        real(8), dimension(3) :: beta, delta
        real(8) :: lambda, chisq_new, tol, reduction_factor, increase_factor
        integer :: i, max_iter, iter
        reduction_factor=0.5
        increase_factor=2
        lambda = 100d0
        tol = 1.0e-14
        max_iter = 1000000

        do iter = 1, max_iter
            ! 计算模型值和残差
            do i = 1, n
                L(i) = quadratic(x(i), A, B, C)
                res(i) = y(i) - L(i)
            end do

            ! 计算雅可比矩阵
            call jacobian(x, y, n, A, B, C, J)

            ! 计算海森矩阵和梯度向量
            H = matmul(transpose(J), J) + lambda * generate_identity_matrix(3)
            beta = matmul(transpose(J), res)

            ! 求解线性方程组 H * delta = beta
            call solve_linear_system(H, beta, delta)

            ! 更新参数
            A = A + delta(1)
            B = B + delta(2)
            C = C + delta(3)

            ! 计算新的误差
            chisq_new = sum(res**2)

            ! 输出当前参数和误差
            !print *, 'Iteration ', iter, ': A = ', A, ', B = ', B, ', C = ', C, ', chisq_new = ', chisq_new

            ! 检查收敛性
            if (abs(chisq_new - chisq_old) < tol) then
                print *, 'Convergence achieved at iteration', iter
                print *, 'chisq_old = ', chisq_old, ' chisq_new = ', chisq_new
                exit
            end if
            ! 调整 lambda 
            if (chisq_new < chisq_old) then 
                lambda = lambda * reduction_factor 
            else 
                lambda = lambda * increase_factor 
             end if
            chisq_old = chisq_new
        end do
    end subroutine lmfit

    subroutine jacobian(x, y, n, A, B, C, J)
        implicit none
        integer, intent(in) :: n
        real(8), dimension(n), intent(in) :: x, y
        real(8), intent(in) :: A, B, C
        real(8), dimension(n, 3), intent(out) :: J
        integer :: i
        real(8) :: L, dL_dA, dL_dB, dL_dC

        do i = 1, n
            L = quadratic(x(i), A, B, C)
            dL_dA = x(i)**2
            dL_dB = x(i)
            dL_dC = 1.0d0
            J(i, 1) = dL_dA
            J(i, 2) = dL_dB
            J(i, 3) = dL_dC
        end do
    end subroutine jacobian

    subroutine solve_linear_system(H, beta, delta)
        implicit none
        real(8), dimension(3, 3), intent(in) :: H
        real(8), dimension(3), intent(in) :: beta
        real(8), dimension(3), intent(out) :: delta
        integer :: ipiv(3), info, n

        n = 3
        delta = beta
        call dgesv(n, 1, H, n, ipiv, delta, n, info)
        if (info /= 0) then
            print *, "Error: LAPACK dgesv failed with error code ", info
        end if
    end subroutine solve_linear_system

    function generate_identity_matrix(n) result(I)
        implicit none
        integer, intent(in) :: n
        integer :: j
        real(8), dimension(n, n) :: I
        I = 0.0d0
        do j = 1, n
            I(j, j) = 1.0d0
        end do
    end function generate_identity_matrix

end program read_and_fit


