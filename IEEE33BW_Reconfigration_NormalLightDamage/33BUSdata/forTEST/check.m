%需要matpower包

%读取标准情况
ieee33=case33bw;
%循环次数9
num_Gen=9999999;
result=zeros(37,num_Gen);
num_Find=0;
for i=0:1:num_Gen
    %生成可通过重构解决的场景，即为连通图
    %随机在case元胞中的branch中11列填入数字，判断是否为连通图，如果是则返回第11列的值
    i=i+1;
    rng_State=round(1*rand(37,1));
    ieee33.branch(:,11)=rng_State;
    ieee33.branch(1,11)=1;
    %检查是否存在孤岛
    check_Result=extract_islands(ieee33);
    if length(check_Result)==1
        num_Find=num_Find+1;
        result(:,num_Find)=ieee33.branch(:,11);
    end
end
%删除重复的列
temp=result';
temp=unique(temp,'rows');
result=temp';

