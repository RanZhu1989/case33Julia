%��Ҫmatpower��

%��ȡ��׼���
ieee33=case33bw;
%ѭ������9
num_Gen=9999999;
result=zeros(37,num_Gen);
num_Find=0;
for i=0:1:num_Gen
    %���ɿ�ͨ���ع�����ĳ�������Ϊ��ͨͼ
    %�����caseԪ���е�branch��11���������֣��ж��Ƿ�Ϊ��ͨͼ��������򷵻ص�11�е�ֵ
    i=i+1;
    rng_State=round(1*rand(37,1));
    ieee33.branch(:,11)=rng_State;
    ieee33.branch(1,11)=1;
    %����Ƿ���ڹµ�
    check_Result=extract_islands(ieee33);
    if length(check_Result)==1
        num_Find=num_Find+1;
        result(:,num_Find)=ieee33.branch(:,11);
    end
end
%ɾ���ظ�����
temp=result';
temp=unique(temp,'rows');
result=temp';

