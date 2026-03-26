robot = robotics.RigidBodyTree;
dhpm=[0 0 1 0;
 0 -pi/2 0.3 0;
 0 0 0.3 0;
 0 pi/2 0 0;
 0 -pi/2 0 0;
 0 0 1.0 0];
bodies = cell(7,1);
joints = cell(7,1);
for i = 1:6
 bodies{i} = rigidBody(['link' num2str(i)]);
 if i==2
 joints{i} = rigidBodyJoint(['jnt' num2str(i)],"prismatic");
 elseif i==3
 joints{i} = rigidBodyJoint(['jnt' num2str(i)],"prismatic");
 else
 joints{i} = rigidBodyJoint(['jnt' num2str(i)],"revolute"); 
 end
 if i==1
 setFixedTransform(joints{i},trvec2tform([0, 0, 0]));
 elseif i==2
 setFixedTransform(joints{i},trvec2tform([0, 0, 1]));
 elseif i==3
 setFixedTransform(joints{i},trvec2tform([0, 0, 1]));
 else
 setFixedTransform(joints{i},trvec2tform([1,0,0]));
 end
 joints{1}.JointAxis=[0,0,1];
 joints{2}.JointAxis=[0,0,1];
 joints{3}.JointAxis=[1,0,0];
 joints{4}.JointAxis=[1,0,0];
 joints{5}.JointAxis=[0,0,1];
 joints{6}.JointAxis=[1,0,0];
 bodies{i}.Joint = joints{i};
 if i == 1 
 addBody(robot,bodies{i},"base")
 else
 addBody(robot,bodies{i},bodies{i-1}.Name)
 end
end
bodies{7}=rigidBody("tool");
joints{7}=rigidBodyJoint("tcp","fixed");
setFixedTransform(joints{7},trvec2tform([0, 0, 0]));
bodies{7}.Joint=joints{7};
addBody(robot,bodies{7},bodies{6}.Name)
showdetails(robot)
show(robot);
figure(Name="Interactive GUI")
gui = interactiveRigidBodyTree(robot,MarkerScaleFactor=0.5);
endEffectorPose = trvec2tform([1, 0.5, 0.5]);
ik = robotics.InverseKinematics('RigidBodyTree', robot);
weights = [1 1 1 0 0 0];
initialGuess = zeros(6,1);
configSol= ik('tool',trvec2tform([1, 0.5, 0.5]) , weights, initialGuess);
disp('Joint angles (radians):');
disp(configSol);
show(robot, configSol);
title('Robot Configuration for IK Solution');